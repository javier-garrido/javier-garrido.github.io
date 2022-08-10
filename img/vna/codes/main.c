/*
 * VNA
 */

#include "ch.h"
#include "hal.h"
#include "usbcfg.h"
#include "si5351.h"
#include "nanovna.h"
#include "fft.h"

#include <chprintf.h>
#include <string.h>
#include <math.h>


static BaseSequentialStream *shell_stream;
bool enable_usart_data = false;

#define ENABLE_TRANSFORM_COMMAND
#define ENABLE_SCANBIN_COMMAND
#define ENABLE_USART_COMMAND

#define VNA_SHELL_FUNCTION(command_name) static void command_name(int argc, char *argv[])

static void apply_CH0_error_term_at(int i);
static void apply_edelay(void);

static uint16_t get_sweep_mask(void);
static void cal_interpolate(void);
static void update_frequencies(bool interpolate);
static int  set_frequency(uint32_t freq);
static void set_frequencies(uint32_t start, uint32_t stop, uint16_t points);
static bool sweep(bool break_on_operation, uint16_t ch_mask);
static void transform_domain(void);

uint8_t sweep_mode = SWEEP_ENABLE;

volatile uint16_t wait_count = 0;
// Punto de barrido actual
static uint16_t p_sweep = 0;
// Buffer I2S
static int16_t rx_buffer[AUDIO_BUFFER_LEN * 2];
// Datos medidos en el barrido
float measured[1][POINTS_COUNT][2];
uint32_t frequencies[POINTS_COUNT];


// Hilo para el sweep
static THD_WORKING_AREA(waThread1, 768);
static THD_FUNCTION(Thread1, arg)
{
  while (1) {
  	wdgReset(&WDGD1);
    bool completed = false;
    if (sweep_mode&(SWEEP_ENABLE|SWEEP_ONCE)) {
      completed = sweep(true, get_sweep_mask());
      sweep_mode&=~SWEEP_ONCE;
    } else {
      __WFI();
    }
    // Process collected data only if scan completed
    if ((sweep_mode & SWEEP_ENABLE) && completed) {
      if (electrical_delay != 0) apply_edelay();
      if ((domain_mode & DOMAIN_MODE) == DOMAIN_TIME) transform_domain();
    }
  }
}


static inline void
pause_sweep(void)
{
  sweep_mode &= ~SWEEP_ENABLE;
}

static inline void
resume_sweep(void)
{
  sweep_mode |= SWEEP_ENABLE;
}

void
toggle_sweep(void)
{
  sweep_mode ^= SWEEP_ENABLE;
}


// Dominio del tiempo
static float
bessel0(float x)
{
  const float eps = 0.0001;

  float ret = 0;
  float term = 1;
  float m = 0;

  while (term  > eps * ret) {
    ret += term;
    ++m;
    term *= (x*x) / (4*m*m);
  }
  return ret;
}

static float
kaiser_window(float k, float n, float beta)
{
  if (beta == 0.0) return 1.0;
  float r = (2 * k) / (n - 1) - 1;
  return bessel0(beta * sqrt(1 - r * r)) / bessel0(beta);
}

static void
transform_domain(void)
{
  float* tmp = (float*)spi_buffer;

  uint16_t window_size = sweep_points, offset = 0;
  uint8_t is_lowpass = FALSE;
  uint8_t td_func = domain_mode & TD_FUNC;
  switch (td_func) {
    case TD_FUNC_BANDPASS:
      offset = 0;
      window_size = sweep_points;
      break;
    case TD_FUNC_LOWPASS_IMPULSE:
    case TD_FUNC_LOWPASS_STEP:
      is_lowpass = TRUE;
      offset = sweep_points;
      window_size = sweep_points * 2;
      break;
  }

  float beta = 0.0f;
  switch (domain_mode & TD_WINDOW) {
    case TD_WINDOW_MINIMUM:
    //beta = 0.0f;
      break;
    case TD_WINDOW_NORMAL:
      beta = 6.0f;
      break;
    case TD_WINDOW_MAXIMUM:
      beta = 13.0f;
      break;
  }

  static float window_scale = 1.0f;
  static uint16_t td_cache = 0;
  uint16_t td_check = (domain_mode & (TD_WINDOW|TD_FUNC))|(sweep_points<<5);
  if (td_cache!=td_check){
    td_cache=td_check;
    if (td_func == TD_FUNC_LOWPASS_STEP)
      window_scale = 1.0f;
    else {
      window_scale = 0.0f;
      for (int i = 0; i < sweep_points; i++)
        window_scale += kaiser_window(i + offset, window_size, beta);
      window_scale = (FFT_SIZE/2) / window_scale;
      if (td_func == TD_FUNC_BANDPASS)
        window_scale *= 2;
    }
  }

  uint16_t ch_mask = get_sweep_mask();
  for (int ch = 0; ch < 1; ch++,ch_mask>>=1) {
    if ((ch_mask&1)==0) continue;
    memcpy(tmp, measured[ch], sizeof(measured[0]));
    for (int i = 0; i < sweep_points; i++) {
      float w = kaiser_window(i + offset, window_size, beta) * window_scale;
      tmp[i * 2 + 0] *= w;
      tmp[i * 2 + 1] *= w;
    }
    for (int i = sweep_points; i < FFT_SIZE; i++) {
      tmp[i * 2 + 0] = 0.0;
      tmp[i * 2 + 1] = 0.0;
    }
    if (is_lowpass) {
      for (int i = 1; i < sweep_points; i++) {
        tmp[(FFT_SIZE - i) * 2 + 0] = tmp[i * 2 + 0];
        tmp[(FFT_SIZE - i) * 2 + 1] = -tmp[i * 2 + 1];
      }
    }

    fft_inverse((float(*)[2])tmp);
    memcpy(measured[ch], tmp, sizeof(measured[0]));
    for (int i = 0; i < sweep_points; i++) {
      measured[ch][i][0] /= (float)FFT_SIZE;
      if (is_lowpass) {
        measured[ch][i][1] = 0.0;
      } else {
        measured[ch][i][1] /= (float)FFT_SIZE;
      }
    }
    if ((domain_mode & TD_FUNC) == TD_FUNC_LOWPASS_STEP) {
      for (int i = 1; i < sweep_points; i++) {
        measured[ch][i][0] += measured[ch][i - 1][0];
      }
    }
  }
}


void set_power(uint8_t value){
  if (value > SI5351_CLK_DRIVE_STRENGTH_8MA) value = SI5351_CLK_DRIVE_STRENGTH_AUTO;
  if (current_props._power == value) return;
  current_props._power = value;
  // Actualizar la potencia
  if (!(sweep_mode&SWEEP_ENABLE)) si5351_set_power(value);
}


static void (*sample_func)(float *gamma) = calculate_gamma;

config_t config = {
  .magic       = CONFIG_MAGIC,
  .dac_value   = 1922,
  ._mode       = VNA_MODE_START_STOP,
  .harmonic_freq_threshold = FREQUENCY_THRESHOLD,
  ._serial_speed = 1000000,
  .bandwidth = BANDWIDTH_1000
};

properties_t current_props;

// Configuración por defecto del VNA
static const trace_t def_trace[TRACES_MAX] = {  //enable, type, channel, reserved, scale, refpos
  { 1, TRC_LOGMAG, 0, 0, 10.0, NGRIDY-1 },
  { 1, TRC_LOGMAG, 1, 0, 10.0, NGRIDY-1 },
  { 1, TRC_SMITH,  0, 0, 1.0, 0 },
  { 1, TRC_PHASE,  1, 0, 90.0, NGRIDY/2 }
};

static const marker_t def_markers[MARKERS_MAX] = {
  { 1, 0, 30, 0 }, { 0, 0, 40, 0 }, { 0, 0, 60, 0 }, { 0, 0, 80, 0 }
};

// Cargar las propiedades por defecto
void load_default_properties(void)
{
  current_props._frequency0   =     50000;    // inicio =  50kHz
  current_props._frequency1   = 900000000;    // fin   = 900MHz
  current_props._sweep_points = POINTS_COUNT_DEFAULT; // Puntos por defecto
  current_props._cal_status   = 0;
  current_props._electrical_delay = 0.0;
  memcpy(current_props._trace, def_trace, sizeof(def_trace));
  memcpy(current_props._markers, def_markers, sizeof(def_markers));
  current_props._velocity_factor = 0.7;
  current_props._active_marker   = 0;
  current_props._domain_mode     = 0;
  current_props._marker_smith_format = MS_RLC;
  current_props._power = SI5351_CLK_DRIVE_STRENGTH_AUTO;
}

int load_properties(uint32_t id)
{
  int r = caldata_recall(id);
  update_frequencies(false);
  return r;
}


// Función de callback I2S DMA
static volatile systime_t ready_time = 0;
void i2s_end_callback(I2SDriver *i2sp, size_t offset, size_t n)
{
  int16_t *p = &rx_buffer[offset];
  (void)i2sp;
  if (wait_count == 0 || chVTGetSystemTimeX() < ready_time) return;
  if (wait_count == config.bandwidth+2)      // Resetear y vaciar buffer
    reset_dsp_accumerator();
  else if (wait_count <= config.bandwidth+1) // Procesamiento de los datos
    dsp_process(p, n);

  --wait_count;
}

static const I2SConfig i2sconfig = {
  NULL,                   // TX Buffer
  rx_buffer,              // RX Buffer
  AUDIO_BUFFER_LEN * 2,   // RX Buffer size
  NULL,                   // tx callback
  i2s_end_callback,       // rx callback
  0,                      // i2scfgr
  0                       // i2spr
};

#define DELAY_CHANNEL_CHANGE   3
#define DELAY_SWEEP_START     50    // Delay para empezar el barrido

#define DSP_START(delay) {ready_time = chVTGetSystemTimeX() + delay; wait_count = config.bandwidth+2;}
#define DSP_WAIT         while (wait_count) {__WFI();}
#define RESET_SWEEP      {p_sweep = 0;}

#define SWEEP_CH0_MEASURE   1

static uint16_t get_sweep_mask(void){
  uint16_t ch_mask = 0;
  int t;
  for (t = 0; t < TRACES_MAX; t++) {
    if (!trace[t].enabled)
      continue;
    if (trace[t].channel == 0) ch_mask|=SWEEP_CH0_MEASURE;
  }
  return ch_mask;
}


// Bucle principal para medir el S11
static bool sweep(bool break_on_operation, uint16_t ch_mask)
{
  int delay;
  if (p_sweep>=sweep_points || break_on_operation == false) RESET_SWEEP;
  if (break_on_operation && ch_mask == 0)
    return false;
  
  // Parpadeo del LED mientras escanea
  palClearPad(GPIOC, GPIOC_LED);

  // Esperar un tiempo para potencia estable
  int st_delay = DELAY_SWEEP_START;
  for (; p_sweep < sweep_points; p_sweep++) {
    delay = set_frequency(frequencies[p_sweep]);
    
    // CH0: REFLEXIÓN, resetear y empezar la medida
    if (ch_mask & SWEEP_CH0_MEASURE) {
      tlv320aic3204_select(0);
      DSP_START(delay+st_delay);
      delay = DELAY_CHANNEL_CHANGE;
      DSP_WAIT;
      (*sample_func)(measured[0][p_sweep]);      // Calcula el coeficiente de reflexión
      if(enable_usart_data) {
		streamWrite((BaseSequentialStream *)&SD1, (uint8_t *)measured[0][p_sweep], sizeof(float)*2);
		shell_printf("%f %f\r\n", measured[0][p_sweep][0], measured[0][p_sweep][1]);
	  }
      if (APPLY_CALIBRATION_AFTER_SWEEP == 0 && (cal_status & CALSTAT_APPLY))
        apply_CH0_error_term_at(p_sweep);
    }

    if (operation_requested && break_on_operation) break;
    st_delay = 0;
  }

  if(enable_usart_data) {
		streamWrite((BaseSequentialStream *)&SD1, (uint8_t *)"\r\n\r\n\r\n\r\n", 8);
		shell_printf("\r\n\r\n\r\n\r\n");
	}

  // Aplica la calibración al final
  if (APPLY_CALIBRATION_AFTER_SWEEP && (cal_status & CALSTAT_APPLY) && p_sweep == sweep_points){
    uint16_t start_sweep;
    for (start_sweep = 0; start_sweep < p_sweep; start_sweep++){
      if (ch_mask & SWEEP_CH0_MEASURE) apply_CH0_error_term_at(start_sweep);
    }
  }

  // Parpadeo del LED mientras escanea
  palSetPad(GPIOC, GPIOC_LED);
  return p_sweep == sweep_points;
}


static int set_frequency(uint32_t freq)
{
  return si5351_set_frequency(freq, current_props._power);
}

void set_bandwidth(uint16_t bw_count){
  config.bandwidth = bw_count&0x1FF;
}

uint32_t get_bandwidth_frequency(uint16_t bw_freq){
  return (AUDIO_ADC_FREQ/AUDIO_SAMPLES_COUNT)/(bw_freq+1);
}

#define MAX_BANDWIDTH      (AUDIO_ADC_FREQ/AUDIO_SAMPLES_COUNT)
#define MIN_BANDWIDTH      ((AUDIO_ADC_FREQ/AUDIO_SAMPLES_COUNT)/512 + 1)


void set_sweep_points(uint16_t points){
  if (points == sweep_points || points > POINTS_COUNT)
    return;

  sweep_points = points;
  update_frequencies(cal_status & CALSTAT_APPLY);

}

#define SCAN_MASK_OUT_FREQ       0b00000001
#define SCAN_MASK_OUT_DATA0      0b00000010
#define SCAN_MASK_OUT_DATA1      0b00000100
#define SCAN_MASK_NO_CALIBRATION 0b00001000
#define SCAN_MASK_BINARY         0b10000000


static void
set_frequencies(uint32_t start, uint32_t stop, uint16_t points)
{
  uint32_t i;
  uint32_t step = (points - 1);
  uint32_t span = stop - start;
  uint32_t delta = span / step;
  uint32_t error = span % step;
  uint32_t f = start, df = step>>1;
  for (i = 0; i <= step; i++, f+=delta) {
    frequencies[i] = f;
    df+=error;
    if (df >=step) {
      f++;
      df -= step;
    }
  }
  // Fuera del rango de barrido
  for (; i < POINTS_COUNT; i++)
    frequencies[i] = 0;
}

static void
update_frequencies(bool interpolate)
{
  uint32_t start, stop;
  start = get_sweep_frequency(ST_START);
  stop  = get_sweep_frequency(ST_STOP);

  set_frequencies(start, stop, sweep_points);
  if (interpolate)
    cal_interpolate();
  RESET_SWEEP;
}

void
set_sweep_frequency(int type, uint32_t freq)
{
  int cal_applied = cal_status & CALSTAT_APPLY;

  // Comprobar si la frecuencia se sale de los límites
  if (type != ST_SPAN && freq < START_MIN)
    freq = START_MIN;
  if (freq > STOP_MAX)
    freq = STOP_MAX;
  uint32_t center, span;
  switch (type) {
    case ST_START:
      config._mode &= ~VNA_MODE_CENTER_SPAN;
      frequency0 = freq;
      if (frequency1 < freq) frequency1 = freq;
      break;
    case ST_STOP:
      config._mode &= ~VNA_MODE_CENTER_SPAN;
      frequency1 = freq;
      if (frequency0 > freq) frequency0 = freq;
      break;
    case ST_CENTER:
      config._mode |= VNA_MODE_CENTER_SPAN;
      center = freq;
      span   = (frequency1 - frequency0)>>1;
      if (span > center - START_MIN)
        span = (center - START_MIN);
      if (span > STOP_MAX - center)
        span = (STOP_MAX - center);
      frequency0 = center - span;
      frequency1 = center + span;
      break;
    case ST_SPAN:
      config._mode |= VNA_MODE_CENTER_SPAN;
      center = (frequency0>>1) + (frequency1>>1);
      span = freq>>1;
      if (center < START_MIN + span)
        center = START_MIN + span;
      if (center > STOP_MAX - span)
        center = STOP_MAX - span;
      frequency0 = center - span;
      frequency1 = center + span;
      break;
    case ST_CW:
      config._mode |= VNA_MODE_CENTER_SPAN;
      frequency0 = freq;
      frequency1 = freq;
      break;
  }
  update_frequencies(cal_applied);
}

uint32_t
get_sweep_frequency(int type)
{
  if (frequency0 > frequency1) {
    uint32_t t = frequency0;
    frequency0 = frequency1;
    frequency1 = t;
  }
  switch (type) {
    case ST_START:  return frequency0;
    case ST_STOP:   return frequency1;
    case ST_CENTER: return frequency0/2 + frequency1/2;
    case ST_SPAN:   return frequency1 - frequency0;
    case ST_CW:     return frequency0;
  }
  return 0;
}


// Corrección de errores

static void
eterm_set(int term, float re, float im)
{
  int i;
  for (i = 0; i < sweep_points; i++) {
    cal_data[term][i][0] = re;
    cal_data[term][i][1] = im;
  }
}

static void
eterm_copy(int dst, int src)
{
  memcpy(cal_data[dst], cal_data[src], sizeof cal_data[dst]);
}


static void
eterm_calc_es(void)
{
  int i;
  for (i = 0; i < sweep_points; i++) {
    // z=1/(jwc*z0) = 1/(2*pi*f*c*z0)
    // s11ao = (z-1)/(z+1) = (1-1/z)/(1+1/z) = (1-jwcz0)/(1+jwcz0)
    // Preparo 1/s11ao para mayor eficiencia
    float c = 50e-15;
    float z0 = 50;
    float z = 2 * VNA_PI * frequencies[i] * c * z0;
    float sq = 1 + z*z;
    float s11aor = (1 - z*z) / sq;
    float s11aoi = 2*z / sq;

    // S11mo'= S11mo - Ed
    // S11ms'= S11ms - Ed
    float s11or = cal_data[CAL_OPEN][i][0] - cal_data[ETERM_ED][i][0];
    float s11oi = cal_data[CAL_OPEN][i][1] - cal_data[ETERM_ED][i][1];
    float s11sr = cal_data[CAL_SHORT][i][0] - cal_data[ETERM_ED][i][0];
    float s11si = cal_data[CAL_SHORT][i][1] - cal_data[ETERM_ED][i][1];
    
    // Es = (S11mo'/s11ao + S11ms')/(S11mo' - S11ms')
    float numr = s11sr + s11or * s11aor - s11oi * s11aoi;
    float numi = s11si + s11oi * s11aor + s11or * s11aoi;
    float denomr = s11or - s11sr;
    float denomi = s11oi - s11si;
    sq = denomr*denomr+denomi*denomi;
    cal_data[ETERM_ES][i][0] = (numr*denomr + numi*denomi)/sq;
    cal_data[ETERM_ES][i][1] = (numi*denomr - numr*denomi)/sq;
  }
  cal_status &= ~CALSTAT_OPEN;
  cal_status |= CALSTAT_ES;
}

static void
eterm_calc_er(int sign)
{
  int i;
  for (i = 0; i < sweep_points; i++) {
    // Er = sign*(1-sign*Es)S11ms'
    float s11sr = cal_data[CAL_SHORT][i][0] - cal_data[ETERM_ED][i][0];
    float s11si = cal_data[CAL_SHORT][i][1] - cal_data[ETERM_ED][i][1];
    float esr = cal_data[ETERM_ES][i][0];
    float esi = cal_data[ETERM_ES][i][1];
    if (sign > 0) {
      esr = -esr;
      esi = -esi;
    }
    esr = 1 + esr;
    float err = esr * s11sr - esi * s11si;
    float eri = esr * s11si + esi * s11sr;
    if (sign < 0) {
      err = -err;
      eri = -eri;
    }
    cal_data[ETERM_ER][i][0] = err;
    cal_data[ETERM_ER][i][1] = eri;
  }
  cal_status &= ~CALSTAT_SHORT;
  cal_status |= CALSTAT_ER;
}

// Et está invertido para mayor eficiencia
static void
eterm_calc_et(void)
{
  int i;
  for (i = 0; i < sweep_points; i++) {
    // Et = 1/(S21mt - Ex)
    float etr = cal_data[CAL_THRU][i][0] - cal_data[CAL_ISOLN][i][0];
    float eti = cal_data[CAL_THRU][i][1] - cal_data[CAL_ISOLN][i][1];
    float sq = etr*etr + eti*eti;
    float invr = etr / sq;
    float invi = -eti / sq;
    cal_data[ETERM_ET][i][0] = invr;
    cal_data[ETERM_ET][i][1] = invi;
  }
  cal_status &= ~CALSTAT_THRU;
  cal_status |= CALSTAT_ET;
}


static void apply_CH0_error_term_at(int i)
{
    // S11m' = S11m - Ed
    // S11a = S11m' / (Er + Es S11m')
    float s11mr = measured[0][i][0] - cal_data[ETERM_ED][i][0];
    float s11mi = measured[0][i][1] - cal_data[ETERM_ED][i][1];
    float err = cal_data[ETERM_ER][i][0] + s11mr * cal_data[ETERM_ES][i][0] - s11mi * cal_data[ETERM_ES][i][1];
    float eri = cal_data[ETERM_ER][i][1] + s11mr * cal_data[ETERM_ES][i][1] + s11mi * cal_data[ETERM_ES][i][0];
    float sq = err*err + eri*eri;
    float s11ar = (s11mr * err + s11mi * eri) / sq;
    float s11ai = (s11mi * err - s11mr * eri) / sq;
    measured[0][i][0] = s11ar;
    measured[0][i][1] = s11ai;
}

static void apply_edelay(void)
{
  int i;
  float real, imag;
  float s, c;
  uint16_t ch_mask = get_sweep_mask();
  for (i=0;i<sweep_points;i++){
    vna_sin_cos(electrical_delay * frequencies[i] * 1E-12, &s, &c);
    if (ch_mask & SWEEP_CH0_MEASURE){
      real = measured[0][i][0];
      imag = measured[0][i][1];
      measured[0][i][0] = real * c - imag * s;
      measured[0][i][1] = imag * c + real * s;
    }
  }
}

void
cal_collect(uint16_t type)
{
  uint16_t dst, src;

  static const struct {
    uint16_t set_flag;
    uint16_t clr_flag;
    uint8_t dst;
    uint8_t src;
 } calibration_set[]={
//    type       set data flag       reset flag              destination source
    [CAL_LOAD] = {CALSTAT_LOAD,  ~(           CALSTAT_APPLY), CAL_LOAD,  0},
    [CAL_OPEN] = {CALSTAT_OPEN,  ~(CALSTAT_ES|CALSTAT_APPLY), CAL_OPEN,  0},
    [CAL_SHORT]= {CALSTAT_SHORT, ~(CALSTAT_ER|CALSTAT_APPLY), CAL_SHORT, 0},
    [CAL_THRU] = {CALSTAT_THRU,  ~(CALSTAT_ET|CALSTAT_APPLY), CAL_THRU,  1},
    [CAL_ISOLN]= {CALSTAT_ISOLN, ~(           CALSTAT_APPLY), CAL_ISOLN, 1},
  };
  if (type >= ARRAY_COUNT(calibration_set)) return;
  cal_status|=calibration_set[type].set_flag;
  cal_status&=calibration_set[type].clr_flag;
  dst = calibration_set[type].dst;
  src = calibration_set[type].src;

  // Hacer un sweep para registrar datos
  uint8_t bw = config.bandwidth;
  if (bw < BANDWIDTH_100)
    config.bandwidth = BANDWIDTH_100;

  sweep(false, SWEEP_CH0_MEASURE);
  config.bandwidth = bw;

  // Copiar datos de calibración
  memcpy(cal_data[dst], measured[src], sizeof measured[0]);
}

void
cal_done(void)
{
  if (!(cal_status & CALSTAT_LOAD))
    eterm_set(ETERM_ED, 0.0, 0.0);
  if ((cal_status & CALSTAT_SHORT) && (cal_status & CALSTAT_OPEN)) {
    eterm_calc_es();
    eterm_calc_er(-1);
  } else if (cal_status & CALSTAT_OPEN) {
    eterm_copy(CAL_SHORT, CAL_OPEN);
    eterm_set(ETERM_ES, 0.0, 0.0);
    eterm_calc_er(1);
  } else if (cal_status & CALSTAT_SHORT) {
    eterm_set(ETERM_ES, 0.0, 0.0);
    cal_status &= ~CALSTAT_SHORT;
    eterm_calc_er(-1);
  } else if (!(cal_status & CALSTAT_ER)){
    eterm_set(ETERM_ER, 1.0, 0.0);
  } else if (!(cal_status & CALSTAT_ES)) {
    eterm_set(ETERM_ES, 0.0, 0.0);
  }
    
  if (!(cal_status & CALSTAT_ISOLN))
    eterm_set(ETERM_EX, 0.0, 0.0);
  if (cal_status & CALSTAT_THRU) {
    eterm_calc_et();
  } else if (!(cal_status & CALSTAT_ET)) {
    eterm_set(ETERM_ET, 1.0, 0.0);
  }

  cal_status |= CALSTAT_APPLY;
}

static void
cal_interpolate(void)
{
  const properties_t *src = caldata_reference();
  uint32_t i, j;
  int eterm;
  if (src == NULL)
    return;

  // Subir no interpolados si hay
  if (frequencies[0] == src->_frequency0 && frequencies[src->_sweep_points-1] == src->_frequency1){
    memcpy(current_props._cal_data, src->_cal_data, sizeof(src->_cal_data));
    cal_status = src->_cal_status;
    return;
  }
  
  uint32_t src_f = src->_frequency0;
  for (i = 0; i < sweep_points; i++) {
    if (frequencies[i] >= src_f)
      break;

    // Rellenar cal_data
    for (eterm = 0; eterm < 5; eterm++) {
      cal_data[eterm][i][0] = src->_cal_data[eterm][0][0];
      cal_data[eterm][i][1] = src->_cal_data[eterm][0][1];
    }
  }

  // Reconstruiur la lista de frecuencias de la fuente
  uint32_t src_points = (src->_sweep_points - 1);
  uint32_t span = src->_frequency1 - src->_frequency0;
  uint32_t delta = span / src_points;
  uint32_t error = span % src_points;
  uint32_t df = src_points>>1;
  j = 0;
  for (; i < sweep_points; i++) {
    uint32_t f = frequencies[i];
    if (f == 0) goto interpolate_finish;
    for (; j < src_points; j++) {
      if (src_f <= f && f < src_f + delta) {
        // Encontrado f entre las frecuencias, en j y j+1
        float k1 = (delta == 0) ? 0.0 : (float)(f - src_f) / delta;
        uint32_t idx = j;
        if (si5351_get_harmonic_lvl(src_f) != si5351_get_harmonic_lvl(src_f+delta)) {
          // f en el armónico anterior, extrapolar de los 2 puntos anteriores
          if (si5351_get_harmonic_lvl(f) == si5351_get_harmonic_lvl(src_f)){
            if (idx >= 1){
              idx--; k1+=1.0;
            }
            else // Límite de puntos
              k1 = 0.0;
          }
          // f en el siguiente armónico, extrapolar de los 2 puntos siguientes
          else {
            if (idx < src_points){
              idx++; k1-=1.0;
            }
            else // Límite de puntos
              k1 = 1.0;
          }
        }
        float k0 = 1.0 - k1;
        for (eterm = 0; eterm < 5; eterm++) {
          cal_data[eterm][i][0] = src->_cal_data[eterm][idx][0] * k0 + src->_cal_data[eterm][idx+1][0] * k1;
          cal_data[eterm][i][1] = src->_cal_data[eterm][idx][1] * k0 + src->_cal_data[eterm][idx+1][1] * k1;
        }
        break;
      }
      df+=error;if (df >=src_points) {src_f++;df -= src_points;}
      src_f+=delta;
    }
    if (j == src_points)
      break;
  }

  // Mayor que la última frecuencia del rango de la fuente
  for (; i < sweep_points; i++) {
    // Rellenar cal_data al final de la fuente
    for (eterm = 0; eterm < 5; eterm++) {
      cal_data[eterm][i][0] = src->_cal_data[eterm][src_points][0];
      cal_data[eterm][i][1] = src->_cal_data[eterm][src_points][1];
    }
  }
interpolate_finish:
  cal_status = src->_cal_status | CALSTAT_INTERPOLATED;
}


void set_electrical_delay(float picoseconds)
{
  if (electrical_delay != picoseconds) {
    electrical_delay = picoseconds;
  }
}


#ifdef ENABLE_TRANSFORM_COMMAND
static void
set_domain_mode(int mode) // DOMAIN_FREQ o DOMAIN_TIME
{
  if (mode != (domain_mode & DOMAIN_MODE)) {
    domain_mode = (domain_mode & ~DOMAIN_MODE) | (mode & DOMAIN_MODE);
    uistat.lever_mode = LM_MARKER;
  }
}

static void
set_timedomain_func(int func) // TD_FUNC_LOWPASS_IMPULSE, TD_FUNC_LOWPASS_STEP o TD_FUNC_BANDPASS
{
  domain_mode = (domain_mode & ~TD_FUNC) | (func & TD_FUNC);
}

static void
set_timedomain_window(int func) // TD_WINDOW_MINIMUM/TD_WINDOW_NORMAL/TD_WINDOW_MAXIMUM
{
  domain_mode = (domain_mode & ~TD_WINDOW) | (func & TD_WINDOW);
}
#endif


/*
 * Shell
 */

#ifdef ENABLE_TRANSFORM_COMMAND
VNA_SHELL_FUNCTION(cmd_transform)
{
  int i;
  if (argc == 0) {
    goto usage;
  }
  //                                         0   1       2    3        4       5      6       7
  static const char cmd_transform_list[] = "on|off|impulse|step|bandpass|minimum|normal|maximum";
  for (i = 0; i < argc; i++) {
    switch (get_str_index(argv[i], cmd_transform_list)) {
      case 0:
        set_domain_mode(DOMAIN_TIME);
        return;
      case 1:
        set_domain_mode(DOMAIN_FREQ);
        return;
      case 2:
        set_timedomain_func(TD_FUNC_LOWPASS_IMPULSE);
        return;
      case 3:
        set_timedomain_func(TD_FUNC_LOWPASS_STEP);
        return;
      case 4:
        set_timedomain_func(TD_FUNC_BANDPASS);
        return;
      case 5:
        set_timedomain_window(TD_WINDOW_MINIMUM);
        return;
      case 6:
        set_timedomain_window(TD_WINDOW_NORMAL);
        return;
      case 7:
        set_timedomain_window(TD_WINDOW_MAXIMUM);
        return;
      default:
        goto usage;
    }
  }
  return;
usage:
  shell_printf("usage: transform {%s} [...]\r\n", cmd_transform_list);
}
#endif

VNA_SHELL_FUNCTION(cmd_freq)
{
  if (argc != 1) {
    goto usage;
  }
  uint32_t freq = my_atoui(argv[0]);

  pause_sweep();
  set_frequency(freq);
  return;
usage:
  shell_printf("usage: freq {frequency(Hz)}\r\n");
}

VNA_SHELL_FUNCTION(cmd_power)
{
  if (argc == 0) {
    shell_printf("power: %d\r\n", current_props._power);
    return;
  }
  if (argc != 1) {
    shell_printf("usage: power {0-3}|{255 - auto}\r\n");
    return;
  }
  set_power(my_atoi(argv[0]));
}

VNA_SHELL_FUNCTION(cmd_threshold)
{
  uint32_t value;
  if (argc != 1) {
    shell_printf("usage: threshold {frequency in harmonic mode}\r\n"\
                 "current: %d\r\n", config.harmonic_freq_threshold);
    return;
  }
  value = my_atoui(argv[0]);
  config.harmonic_freq_threshold = value;
}

VNA_SHELL_FUNCTION(cmd_saveconfig)
{
  (void)argc;
  (void)argv;
  config_save();
  shell_printf("Config saved.\r\n");
}

VNA_SHELL_FUNCTION(cmd_clearconfig)
{
  if (argc != 1) {
    shell_printf("usage: clearconfig {protection key}\r\n");
    return;
  }

  if (strcmp(argv[0], "1234") != 0) {
    shell_printf("Key unmatched.\r\n");
    return;
  }

  clear_all_config_prop_data();
  shell_printf("Config and all cal data cleared.\r\n"\
               "Do reset manually to take effect.\r\n");
}

VNA_SHELL_FUNCTION(cmd_pause)
{
  (void)argc;
  (void)argv;
  pause_sweep();
}

VNA_SHELL_FUNCTION(cmd_resume)
{
  (void)argc;
  (void)argv;
  update_frequencies(cal_status & CALSTAT_APPLY);
  resume_sweep();
}

VNA_SHELL_FUNCTION(cmd_reset)
{
  (void)argc;
  (void)argv;
  if (argc == 1) {
    if (strcmp(argv[0], "dfu") == 0) {
      shell_printf("Performing reset to DFU mode\r\n");
      enter_dfu();
      return;
    }
  }
  shell_printf("Performing reset\r\n");
  rccEnableWWDG(FALSE);
  WWDG->CFR = 0x60;
  WWDG->CR = 0xff;
  while (1);
}

VNA_SHELL_FUNCTION(cmd_data)
{
  int i;
  int sel = 0;
  float (*array)[2];
  if (argc == 1)
    sel = my_atoi(argv[0]);
  if (sel < 0 || sel >=7)
    goto usage;

  array = sel < 2 ? measured[sel] : cal_data[sel-2];

  for (i = 0; i < sweep_points; i++)
    shell_printf("%f %f\r\n", array[i][0], array[i][1]);
  return;
usage:
  shell_printf("usage: data [array]\r\n");
}

VNA_SHELL_FUNCTION(cmd_bandwidth)
{
  uint16_t user_bw;
  if (argc == 1)
    user_bw = my_atoui(argv[0]);
  else if (argc == 2){
    uint16_t f = my_atoui(argv[0]);
         if (f > MAX_BANDWIDTH) user_bw = 0;
    else if (f < MIN_BANDWIDTH) user_bw = 511;
    else user_bw = ((AUDIO_ADC_FREQ+AUDIO_SAMPLES_COUNT/2)/AUDIO_SAMPLES_COUNT)/f - 1;
  }
  else
    goto result;
  set_bandwidth(user_bw);
result:
  shell_printf("bandwidth %d (%uHz)\r\n", config.bandwidth, get_bandwidth_frequency(config.bandwidth));
}

VNA_SHELL_FUNCTION(cmd_scan)
{
  uint32_t start, stop;
  uint16_t points = sweep_points;
  int i;
  if (argc < 2 || argc > 4) {
    shell_printf("usage: scan {start(Hz)} {stop(Hz)} [points] [outmask]\r\n");
    return;
  }

  start = my_atoui(argv[0]);
  stop = my_atoui(argv[1]);
  if (start == 0 || stop == 0 || start > stop) {
      shell_printf("frequency range is invalid\r\n");
      return;
  }
  if (argc >= 3) {
    points = my_atoui(argv[2]);
    if (points == 0 || points > POINTS_COUNT) {
      shell_printf("sweep points exceeds range "define_to_STR(POINTS_COUNT)"\r\n");
      return;
    }
    sweep_points = points;
  }
  uint16_t mask = 0;
  uint16_t sweep_ch = SWEEP_CH0_MEASURE;

#ifdef ENABLE_SCANBIN_COMMAND
  if (argc == 4) {
    mask = my_atoui(argv[3]);
    if (sweep_mode&SWEEP_BINARY) mask|=SCAN_MASK_BINARY;
    sweep_ch = (mask>>1)&3;
  }
  sweep_mode&=~(SWEEP_BINARY);
#endif

  uint32_t old_cal_status = cal_status;
  if (mask&SCAN_MASK_NO_CALIBRATION) cal_status&=~CALSTAT_APPLY;
  // Reconstruir la tabla de frecuencias si es necesario
  if (frequencies[0]!=start || frequencies[points-1]!=stop){
    set_frequencies(start, stop, points);
    if (cal_status & CALSTAT_APPLY)
      cal_interpolate();
  }

  if (sweep_ch & SWEEP_CH0_MEASURE)
    sweep(false, sweep_ch);

  cal_status = old_cal_status;

  pause_sweep();
  // Datos de salida
  if (mask) {
    if (mask&SCAN_MASK_BINARY){
      streamWrite(shell_stream, (void *)&mask, sizeof(uint16_t));
      streamWrite(shell_stream, (void *)&points, sizeof(uint16_t));
      for (i = 0; i < points; i++) {
        if (mask & SCAN_MASK_OUT_FREQ ) streamWrite(shell_stream, (void *)&frequencies[i],    sizeof(uint32_t));  // 4 bytes .. frecuencia
        if (mask & SCAN_MASK_OUT_DATA0) streamWrite(shell_stream, (void *)&measured[0][i][0], sizeof(float)* 2);  // 4+4 bytes .. S11 real/imag
      }
    }
    else{
      for (i = 0; i < points; i++) {
        if (mask & SCAN_MASK_OUT_FREQ ) shell_printf("%u ", frequencies[i]);
        if (mask & SCAN_MASK_OUT_DATA0) shell_printf("%f %f ", measured[0][i][0], measured[0][i][1]);
        shell_printf("\r\n");
      }
    }
  }
}

#ifdef ENABLE_SCANBIN_COMMAND
VNA_SHELL_FUNCTION(cmd_scan_bin)
{
  sweep_mode|= SWEEP_BINARY;
  cmd_scan(argc, argv);
  sweep_mode&=~(SWEEP_BINARY);
}
#endif

VNA_SHELL_FUNCTION(cmd_sweep)
{
  if (argc == 0) {
    shell_printf("%u %u %d\r\n", get_sweep_frequency(ST_START), get_sweep_frequency(ST_STOP), sweep_points);
    return;
  } else if (argc > 3) {
    goto usage;
  }
  uint32_t value0 = 0;
  uint32_t value1 = 0;
  uint32_t value2 = 0;
  if (argc >= 1) value0 = my_atoui(argv[0]);
  if (argc >= 2) value1 = my_atoui(argv[1]);
  if (argc >= 3) value2 = my_atoui(argv[2]);
#if MAX_FREQ_TYPE != 5
#error "Sweep mode possibly changed, check cmd_sweep function"
#endif
  // Comprobar sweep {start|stop|center|span|cw} {freq(Hz)}
  // Obtengo ST_START, ST_STOP, ST_CENTER, ST_SPAN, ST_CW
  static const char sweep_cmd[] = "start|stop|center|span|cw";
  if (argc == 2 && value0 == 0) {
    int type = get_str_index(argv[0], sweep_cmd);
    if (type == -1)
      goto usage;
    set_sweep_frequency(type, value1);
    return;
  }
  //  Comprobar sweep {start(Hz)} [stop(Hz)]
  if (value0)
    set_sweep_frequency(ST_START, value0);
  if (value1)
    set_sweep_frequency(ST_STOP, value1);
  if (value2)
    set_sweep_points(value2);
  return;
usage:
  shell_printf("usage: sweep {start(Hz)} [stop(Hz)] [points]\r\n"\
               "\tsweep {%s} {freq(Hz)}\r\n", sweep_cmd);
}

VNA_SHELL_FUNCTION(cmd_cal)
{
  static const char *items[] = { "load", "open", "short", "thru", "isoln", "Es", "Er", "Et", "cal'ed" };

  if (argc == 0) {
    int i;
    for (i = 0; i < 9; i++) {
      if (cal_status & (1<<i))
        shell_printf("%s ", items[i]);
    }
    shell_printf("\r\n");
    return;
  }
  //                                     0    1     2    3     4    5  6   7     8
  static const char cmd_cal_list[] = "load|open|short|thru|isoln|done|on|off|reset";
  switch (get_str_index(argv[0], cmd_cal_list)) {
    case 0:
      cal_collect(CAL_LOAD);
      return;
    case 1:
      cal_collect(CAL_OPEN);
      return;
    case 2:
      cal_collect(CAL_SHORT);
      return;
    case 3:
      cal_collect(CAL_THRU);
      return;
    case 4:
      cal_collect(CAL_ISOLN);
      return;
    case 5:
      cal_done();
      return;
    case 6:
      cal_status |= CALSTAT_APPLY;
      return;
    case 7:
      cal_status &= ~CALSTAT_APPLY;
      return;
    case 8:
      cal_status = 0;
      return;
    default:
      break;
  }
  shell_printf("usage: cal [%s]\r\n", cmd_cal_list);
}

VNA_SHELL_FUNCTION(cmd_save)
{
  if (argc != 1)
    goto usage;

  int id = my_atoi(argv[0]);
  if (id < 0 || id >= SAVEAREA_MAX)
    goto usage;
  caldata_save(id);
  return;

 usage:
  shell_printf("save {id}\r\n");
}

VNA_SHELL_FUNCTION(cmd_recall)
{
  if (argc != 1)
    goto usage;

  int id = my_atoi(argv[0]);
  if (id < 0 || id >= SAVEAREA_MAX)
    goto usage;
  if (load_properties(id))
    shell_printf("Err, default load\r\n");
  return;
 usage:
  shell_printf("recall {id}\r\n");
}

VNA_SHELL_FUNCTION(cmd_edelay)
{
  if (argc != 1) {
    shell_printf("%f\r\n", electrical_delay);
    return;
  }
  set_electrical_delay(my_atof(argv[0]));
}

VNA_SHELL_FUNCTION(cmd_frequencies)
{
  int i;
  (void)argc;
  (void)argv;
  for (i = 0; i < sweep_points; i++) {
    if (frequencies[i] == 0) break;
    shell_printf("%u\r\n", frequencies[i]);
  }
}

#ifdef __USE_SERIAL_CONSOLE__
#ifdef ENABLE_USART_COMMAND
VNA_SHELL_FUNCTION(cmd_usart_cfg)
{
  if (argc != 1) goto result;
  uint32_t speed = my_atoui(argv[0]);
  if (speed < 300) speed = 300;
  config._serial_speed = speed;
  shell_update_speed();
result:
  shell_printf("Serial: %u baud\r\n", config._serial_speed);
}

VNA_SHELL_FUNCTION(cmd_usart)
{
  uint32_t time = 2000; // 200ms de espera por defecto
  if (argc == 0 || argc > 2 || (config._mode & VNA_MODE_SERIAL)) return;
  if (argc == 2) time = my_atoui(argv[1])*10;
  sdWriteTimeout(&SD1, (uint8_t *)argv[0], strlen(argv[0]), time);
  sdWriteTimeout(&SD1, (uint8_t *)VNA_SHELL_NEWLINE_STR, sizeof(VNA_SHELL_NEWLINE_STR)-1, time);
  uint32_t size;
  uint8_t buffer[64];
  while ((size = sdReadTimeout(&SD1, buffer, sizeof(buffer), time)))
    streamWrite(&SDU1, buffer, size);
}

VNA_SHELL_FUNCTION(cmd_usart_data)
{
  if(argc == 1) {
  uint8_t resp = my_atoui(argv[0]);
    if(resp == 1) {
      enable_usart_data = true;
      return;
    }
    if(resp == 0) {
      enable_usart_data = false;
      return;
    }
  }
  shell_printf("usage: usart_data {1/0}\r\n");
}
#endif
#endif

VNA_SHELL_FUNCTION(cmd_help)
{
  (void)argc;
  (void)argv;
  const VNAShellCommand *scp = commands;
  shell_printf("Commands:");
  while (scp->sc_name != NULL) {
    shell_printf(" %s", scp->sc_name);
    scp++;
  }
  shell_printf(VNA_SHELL_NEWLINE_STR);
  return;
}


#pragma pack(push, 2)
typedef struct {
  const char           *sc_name;
  vna_shellcmd_t    sc_function;
  uint16_t flags;
} VNAShellCommand;
#pragma pack(pop)

// Algunos comandos solo se pueden ejecutar en el hilo del sweep, no en el main
#define CMD_WAIT_MUTEX  1
#define CMD_BREAK_SWEEP 2

// Lista de comandos
static const VNAShellCommand commands[] =
{
    {"scan"        , cmd_scan        , CMD_WAIT_MUTEX|CMD_BREAK_SWEEP},
#ifdef ENABLE_SCANBIN_COMMAND
    {"scan_bin"    , cmd_scan_bin    , CMD_WAIT_MUTEX|CMD_BREAK_SWEEP},
#endif
    {"data"        , cmd_data        , 0},
    {"frequencies" , cmd_frequencies , 0},
    {"freq"        , cmd_freq        , CMD_WAIT_MUTEX|CMD_BREAK_SWEEP},
    {"sweep"       , cmd_sweep       , CMD_WAIT_MUTEX|CMD_BREAK_SWEEP},
    {"power"       , cmd_power       , 0},
    {"bandwidth"   , cmd_bandwidth   , 0},
    {"saveconfig"  , cmd_saveconfig  , 0},
    {"clearconfig" , cmd_clearconfig , 0},
    {"pause"       , cmd_pause       , CMD_WAIT_MUTEX|CMD_BREAK_SWEEP},
    {"resume"      , cmd_resume      , CMD_WAIT_MUTEX|CMD_BREAK_SWEEP},
    {"cal"         , cmd_cal         , CMD_WAIT_MUTEX},
    {"save"        , cmd_save        , 0},
    {"recall"      , cmd_recall      , CMD_WAIT_MUTEX|CMD_BREAK_SWEEP},
    {"edelay"      , cmd_edelay      , 0},
    {"reset"       , cmd_reset       , 0},
#ifdef __USE_SERIAL_CONSOLE__
#ifdef ENABLE_USART_COMMAND
    {"usart_cfg"   , cmd_usart_cfg   , CMD_WAIT_MUTEX|CMD_BREAK_SWEEP},
    {"usart"       , cmd_usart       , CMD_WAIT_MUTEX|CMD_BREAK_SWEEP},
    {"usart_data"  , cmd_usart_data  , 0},
#endif
#endif
#ifdef ENABLE_TRANSFORM_COMMAND
    {"transform"   , cmd_transform   , 0},
#endif
    {"threshold"   , cmd_threshold   , 0},
    {"help"        , cmd_help        , 0},
    {NULL          , NULL            , 0}
};


// Comprobar conexión Serial
#ifdef __USE_SERIAL_CONSOLE__
#if HAL_USE_SERIAL == FALSE
#error "For serial console need HAL_USE_SERIAL as TRUE in halconf.h"
#endif

// Seleccionar stream de entrada para comandos de la Shell
#define PREPARE_STREAM shell_stream = (config._mode&VNA_MODE_SERIAL) ? (BaseSequentialStream *)&SD1 : (BaseSequentialStream *)&SDU1;

// Actualizar la velocidad y configuración de la conexión Serial
void shell_update_speed(void){
  SerialConfig s_config = {config._serial_speed, 0, USART_CR2_STOP1_BITS, 0 };
  sdStop(&SD1);
  sdStart(&SD1, &s_config);  // USART
}

// Comprobar el estado de la conexión USB
static bool usb_IsActive(void){
  return usbGetDriverStateI(&USBD1) == USB_ACTIVE;
}

// Resetear la cola I/O de la Shell
void shell_reset_console(void){
  // USB
  if (usb_IsActive()){
    if (config._mode & VNA_MODE_SERIAL)
      sduDisconnectI(&SDU1);
    else
      sduConfigureHookI(&SDU1);
  }
  // Serial
  oqResetI(&SD1.oqueue);
  iqResetI(&SD1.iqueue);
}

// Comprobar conexión activa para la Shell
static bool shell_check_connect(void){
  // Serial
  if (config._mode & VNA_MODE_SERIAL)
    return true;
  // USB
  return usb_IsActive();
}

static void shell_init_connection(void){
  // Inicializa y comienza el driver serial-over-USB CDC SDU1, conectado a USBD1
  sduObjectInit(&SDU1);
  sduStart(&SDU1, &serusbcfg);
  
  // Establecer la configuración de la velocidad Serial para SD1
  shell_update_speed();

  // Activa el driver del USB y luego el pull-up del bus del USB en D+
  usbDisconnectBus(&USBD1);
  chThdSleepMilliseconds(100);
  usbStart(&USBD1, &usbcfg);
  usbConnectBus(&USBD1);

  // Establecer el stream I/O (SDU1 o SD1) para la Shell
  PREPARE_STREAM;
}

#else

// Solo consola por USB, shell_stream siempre sobre el USB
#define PREPARE_STREAM

// Comprobar si la conexión está activa
static bool shell_check_connect(void){
  return SDU1.config->usbp->state == USB_ACTIVE;
}

// Inicializar la conexión I/O con la Shell sobre el USB
static void shell_init_connection(void){
  // Inicializa y comienza el driver serial-over-USB CDC SDU1, conectado a USBD1
  sduObjectInit(&SDU1);
  sduStart(&SDU1, &serusbcfg);

  // Activa el driver del USB y luego el pull-up del bus del USB en D+
  usbDisconnectBus(&USBD1);
  chThdSleepMilliseconds(100);
  usbStart(&USBD1, &usbcfg);
  usbConnectBus(&USBD1);

  // Establecer el stream I/O SDU1 para la Shell
  shell_stream = (BaseSequentialStream *)&SDU1;
}
#endif


/*
 * Leer comando de shell_stream
 */
static int VNAShell_readLine(char *line, int max_size)
{
  // Leer la línea del stream de entrada
  uint8_t c;
  // Preparar I/O para shell_stream
  PREPARE_STREAM;
  char *ptr = line;
  while (1) {
    // Devuelve 0 solo si el stream no está activo
    if (streamRead((BaseSequentialStream *)&SD1, &c, 1) == 0)
      return 0;
    // Borrar
    if (c == 8 || c == 0x7f) {
      if (ptr != line) {
        static const char backspace[] = {0x08, 0x20, 0x08, 0x00};
        shell_printf(backspace);
        ptr--;
      }
      continue;
    }
    // Nueva línea (Enter)
    if (c == '\r') {
      shell_printf(VNA_SHELL_NEWLINE_STR);
      *ptr = 0;
      return 1;
    }
    // Otros (los salta)
    if (c < 0x20)
      continue;
    // Almacenar
    if (ptr < line + max_size - 1) {
      streamPut(shell_stream, c); 		// Echo
      *ptr++ = (char)c;
    }
  }
  return 0;
}

/*
 * Comprobar y ejecutar el comando
 */
static void VNAShell_executeLine(char *line)
{
  // Comprueba y ejecuta la línea
  char *lp = line, *ep;
  shell_nargs = 0;

  while (*lp != 0) {
    // Se salta los espacios en blanco y tabulaciones al comienzo del string
    while (*lp == ' ' || *lp == '\t') lp++;
    // Si un argumento comienza con comillas, entonces su delimitador son otras comillas,
    // si no, el delimitador es un espacio en blanco
    ep = (*lp == '"') ? strpbrk(++lp, "\"") : strpbrk(lp, " \t");
    // Almacenar en la string de args
    shell_args[shell_nargs++] = lp;
    // Final de la string de entrada
    if ((lp = ep) == NULL) break;
    // Comprobación de los límites del argumento
    if (shell_nargs > VNA_SHELL_MAX_ARGUMENTS) {
      shell_printf("too many arguments, max " define_to_STR(
          VNA_SHELL_MAX_ARGUMENTS) "" VNA_SHELL_NEWLINE_STR);
      return;
    }
    // Establecer cero al final de la string y continuar comprobando
    *lp++ = 0;
  }
  if (shell_nargs == 0) return;
  // Ejecutar línea
  const VNAShellCommand *scp;
  for (scp = commands; scp->sc_name != NULL; scp++) {
    if (strcmp(scp->sc_name, shell_args[0]) == 0) {
      if (scp->flags & CMD_WAIT_MUTEX) {
        shell_function = scp->sc_function;
        if (scp->flags & CMD_BREAK_SWEEP) operation_requested|=OP_CONSOLE;
        // Espera a ejecutar el comando en el hilo del sweep
        do {
          chThdSleepMilliseconds(100);
        } while (shell_function);
      } else {
        scp->sc_function(shell_nargs - 1, &shell_args[1]);
      }
      return;
    }
  }
  shell_printf("%s?" VNA_SHELL_NEWLINE_STR, shell_args[0]);
}


/*
 * Configuración del bus I2C
 */

// Define la velocidad del bus I2C
#define STM32_I2C_SPEED                     600

#if STM32_I2C1_CLOCK == 8    // STM32_I2C1SW == STM32_I2C1SW_HSI     (HSI=8MHz)
#if   STM32_I2C_SPEED == 400 // 400kHz @ HSI 8MHz
 #define STM32_I2C_TIMINGR  STM32_TIMINGR_PRESC(0U)  |\
                            STM32_TIMINGR_SCLDEL(3U) | STM32_TIMINGR_SDADEL(1U) |\
                            STM32_TIMINGR_SCLH(3U)   | STM32_TIMINGR_SCLL(9U)
#endif
#elif  STM32_I2C1_CLOCK == 48 // STM32_I2C1SW == STM32_I2C1SW_SYSCLK  (SYSCLK = 48MHz)
 #if   STM32_I2C_SPEED == 400 // 400kHz @ SYSCLK 48MHz
 #define STM32_I2C_TIMINGR  STM32_TIMINGR_PRESC(5U)  |\
                            STM32_TIMINGR_SCLDEL(3U) | STM32_TIMINGR_SDADEL(3U) |\
                            STM32_TIMINGR_SCLH(3U)   | STM32_TIMINGR_SCLL(9U)
 #elif STM32_I2C_SPEED == 600 // 600kHz @ SYSCLK 48MHz, obtener valores manualmente, x1.5 velocidad I2C
 #define STM32_I2C_TIMINGR  STM32_TIMINGR_PRESC(0U)   |\
                            STM32_TIMINGR_SCLDEL(10U) | STM32_TIMINGR_SDADEL(10U) |\
                            STM32_TIMINGR_SCLH(30U)   | STM32_TIMINGR_SCLL(50U)
 #elif STM32_I2C_SPEED == 900 // 900kHz @ SYSCLK 48MHz, obtener valores manualmente, x2 velocidad I2C
 #define STM32_I2C_TIMINGR  STM32_TIMINGR_PRESC(0U)   |\
                            STM32_TIMINGR_SCLDEL(10U) | STM32_TIMINGR_SDADEL(10U) |\
                            STM32_TIMINGR_SCLH(23U)   | STM32_TIMINGR_SCLL(30U)
 #endif
#elif  STM32_I2C1_CLOCK == 72 // STM32_I2C1SW == STM32_I2C1SW_SYSCLK  (SYSCLK = 72MHz)
 #if   STM32_I2C_SPEED == 400 // ~400kHz @ SYSCLK 72MHz
 #define STM32_I2C_TIMINGR  STM32_TIMINGR_PRESC(0U)   |\
                            STM32_TIMINGR_SCLDEL(10U) | STM32_TIMINGR_SDADEL(10U) |\
                            STM32_TIMINGR_SCLH(80U)   | STM32_TIMINGR_SCLL(100U)
 #elif STM32_I2C_SPEED == 600 // ~600kHz @ SYSCLK 72MHz, obtener valores manualmente, x1.5 velocidad I2C
 #define STM32_I2C_TIMINGR  STM32_TIMINGR_PRESC(0U)   |\
                            STM32_TIMINGR_SCLDEL(10U) | STM32_TIMINGR_SDADEL(10U) |\
                            STM32_TIMINGR_SCLH(40U)   | STM32_TIMINGR_SCLL(80U)
 #elif STM32_I2C_SPEED == 900 // ~900kHz @ SYSCLK 72MHz, obtener valores manualmente, x2 velocidad I2C
 #define STM32_I2C_TIMINGR  STM32_TIMINGR_PRESC(0U)   |\
                            STM32_TIMINGR_SCLDEL(10U) | STM32_TIMINGR_SDADEL(10U) |\
                            STM32_TIMINGR_SCLH(30U)   | STM32_TIMINGR_SCLL(40U)
 #endif
#endif


// Configuración del reloj I2C (depende de STM32_I2C1SW en mcuconf.h)
static const I2CConfig i2ccfg = {
  .timingr = STM32_I2C_TIMINGR,  // Inicialización del registro TIMINGR
  .cr1 = 0,                      // Inicialización del registro CR1
  .cr2 = 0                       // Inicialización del registro CR2
};


/*
 * Watchdog frecuencia (más de un segundo): LSI = 40000 / (64 * 1000)
 */
static const WDGConfig wdgcfg = {
	STM32_IWDG_PR_64,
	STM32_IWDG_RL(1250),
	STM32_IWDG_WIN_DISABLED
};


/*
 * Main
 */
// El tamaño del stack del hilo del main está definido en el Makefile, USE_PROCESS_STACKSIZE = 0x200
int main(void)
{
/*
 * Inicializar el sistema ChibiOS
 */
  halInit();
  chSysInit();

/*
 * Restaurar la configuración
 */
  config_recall();

/*
 * Restaurar las frecuencias y la calibración del slot 0 de la memoria flash
 */
  load_properties(0);

/*
 * Inicializar la conexión con la Shell
 */
  shell_init_connection();

/*
 * Inicializar el I2C
 */
  i2cStart(&I2CD1, &i2ccfg);

/*
 * Inicializar el Si5351
 */
  si5351_init();

/*
 * Inicializar el códec TLV320AIC3204
 */
  tlv320aic3204_init();

/*
 * Inicializar el I2S
 */
  i2sInit();
  i2sObjectInit(&I2SD2);
  i2sStart(&I2SD2, &i2sconfig);
  i2sStartExchange(&I2SD2);

 /*
 * Inicializar el Watchdog
 */
  wdgStart(&WDGD1, &wdgcfg);

/*
 * Comienzo del hilo sweep
 */
  chThdCreateStatic(waThread1, sizeof(waThread1), NORMALPRIO-1, Thread1, NULL);

/*
 * Lectura y ejecución de comandos de la Shell
 */
  while (1) {
    if (shell_check_connect()) {
      shell_printf(VNA_SHELL_NEWLINE_STR"NanoVNA Shell"VNA_SHELL_NEWLINE_STR);

      do {
        shell_printf(VNA_SHELL_PROMPT_STR);
        if (VNAShell_readLine(shell_line, VNA_SHELL_MAX_LENGTH))
          VNAShell_executeLine(shell_line);
        else
          chThdSleepMilliseconds(200);
      } while (shell_check_connect());

    }
    chThdSleepMilliseconds(1000);
  }
}