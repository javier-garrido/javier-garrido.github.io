from machine import Pin, UART
import time
import struct
import math
import cmath
import os

DEBUG = False

# Comandos para guardar un barrido por UART
mensajeON = 'usart_data 1\r\n'
mensajeOFF = 'usart_data 0\r\n'

# Comandos para calibrar
mensajeCALL = 'cal load\r\n'
mensajeCALO = 'cal open\r\n'
mensajeCALS = 'cal short\r\n'

# Dirección del VNA para SDI-12
ID = '1'

# UART para la comunicación con el VNA
uart = UART(2, 1000000)
uart.init(1000000, bits=8, parity=None, stop=1)

# UART para SDI-12
uart_sdi12 = UART(1, 1200)
uart_sdi12.init(1200, bits=7, parity=0, stop=1)
buffer_sdi12 = b''

botonON = Pin(12, Pin.IN)
botonCALL = Pin(13, Pin.IN)
botonCALO = Pin(14, Pin.IN)
botonCALS = Pin(27, Pin.IN)
direccion_sdi12 = Pin(15, Pin.OUT)

# Comenzar barrido
def pulsarON(p):
    for i in range(len(mensajeON)):
        if DEBUG:
            print(mensajeON[i], end='')
        uart.write(str(mensajeON[i]))
        if DEBUG:
            time.sleep(0.1)

# Calibración "load"
def pulsarCALL(p):
    for i in range(len(mensajeCALL)):
        if DEBUG:
            print(mensajeCALL[i], end='')
        uart.write(str(mensajeCALL[i]))
        if DEBUG:
            time.sleep(0.1)

# Calibración "open"
def pulsarCALO(p):
    for i in range(len(mensajeCALO)):
        if DEBUG:
            print(mensajeCALO[i], end='')
        uart.write(str(mensajeCALO[i]))
        if DEBUG:
            time.sleep(0.1)

# Calibración "short"
def pulsarCALS(p):
    for i in range(len(mensajeCALS)):
        if DEBUG:
            print(mensajeCALS[i], end='')
        uart.write(str(mensajeCALS[i]))
        if DEBUG:
            time.sleep(0.1)

botonON.irq(trigger=Pin.IRQ_RISING, handler=pulsarON)
botonCALL.irq(trigger=Pin.IRQ_RISING, handler=pulsarCALL)
botonCALO.irq(trigger=Pin.IRQ_RISING, handler=pulsarCALO)
botonCALS.irq(trigger=Pin.IRQ_RISING, handler=pulsarCALS)

measured = [[None for i in range(2)] for j in range(101)]
permitividad_compleja = [None for i in range(101)]
permitividad_re_im = [[None for i in range(2)] for j in range(101)]
i = 0
medido = 0
num_fichero = '1'
n_puntos = 5

# Constantes para la permitividad
c =  [[complex(-1.69, -7.32), complex(-1.79, -7.49), complex(1.00, -0.01)], [complex(-83.93, -867.51), complex(-137.77, -860.57), complex(-0.85, 0.23)],
      [complex(13.16, -460.70), complex(-44.05, -458.85), complex(-1.00, 0.07)], [complex(17.03, -313.73), complex(-41.43, -311.76), complex(-1.06, 0.11)],
      [complex(2.09, -234.08), complex(-55.36, -228.01), complex(-0.99, 0.25)], [complex(-1.70, -183.19), complex(-57.20, -174.37), complex(-0.92, 0.32)],
      [complex(0.92, -153.68), complex(-55.50, -143.64), complex(-0.95, 0.37)], [complex(1.61, -131.26), complex(-54.23, -120.17), complex(-0.93, 0.42)],
      [complex(0.85, -114.67), complex(-54.09, -101.54), complex(-0.89, 0.48)], [complex(0.06, -102.54), complex(-54.40, -87.25), complex(-0.85, 0.55)],
      [complex(-0.13, -92.08), complex(-53.80, -75.13), complex(-0.81, 0.60)], [complex(-0.51, -83.48), complex(-53.21, -64.79), complex(-0.76, 0.66)],
      [complex(-0.61, -75.83), complex(-51.74, -55.61), complex(-0.71, 0.69)], [complex(-0.49, -70.63), complex(-51.14, -48.77), complex(-0.67, 0.75)],
      [complex(-0.73, -65.68), complex(-50.50, -42.19), complex(-0.62, 0.79)], [complex(-1.14, -60.98), complex(-49.41, -35.69), complex(-0.56, 0.82)],
      [complex(-0.84, -57.32), complex(-48.42, -30.81), complex(-0.51, 0.86)], [complex(-1.03, -54.11), complex(-47.55, -25.94), complex(-0.45, 0.89)],
      [complex(-1.19, -50.87), complex(-46.21, -21.34), complex(-0.39, 0.91)], [complex(-1.00, -48.25), complex(-45.01, -17.70), complex(-0.33, 0.94)],
      [complex(0.14, -46.72), complex(-43.34, -16.10), complex(-0.29, 0.97)], [complex(-0.93, -43.56), complex(-41.54, -11.44), complex(-0.21, 0.94)],
      [complex(-0.85, -41.59), complex(-40.22, -8.71), complex(-0.15, 0.97)], [complex(-1.94, -39.46), complex(-38.69, -4.92), complex(-0.07, 0.95)],
      [complex(-1.59, -38.02), complex(-37.23, -2.71), complex(-0.02, 0.97)], [complex(-1.67, -36.22), complex(-35.54, -0.26), complex(0.04, 0.96)],
      [complex(2.04, -32.01), complex(-31.02, -2.89), complex(-0.04, 1.01)], [complex(-2.20, -33.34), complex(-32.42, 3.99), complex(0.16, 0.95)],
      [complex(-2.62, -32.20), complex(-30.78, 5.94), complex(0.22, 0.93)], [complex(-2.59, -31.35), complex(-29.64, 7.61), complex(0.28, 0.91)],
      [complex(-2.29, -30.20), complex(-28.12, 8.64), complex(0.33, 0.92)], [complex(-2.57, -28.82), complex(-26.48, 10.12), complex(0.40, 0.89)],
      [complex(-3.05, -28.21), complex(-25.14, 11.76), complex(0.47, 0.86)], [complex(-3.00, -27.01), complex(-23.50, 12.70), complex(0.52, 0.82)],
      [complex(-3.13, -25.98), complex(-21.75, 13.85), complex(0.57, 0.78)], [complex(-3.40, -24.52), complex(-19.46, 14.57), complex(0.61, 0.72)],
      [complex(-3.42, -23.66), complex(-17.89, 15.33), complex(0.66, 0.68)], [complex(-3.49, -22.83), complex(-16.27, 15.88), complex(0.71, 0.62)],
      [complex(-3.43, -22.19), complex(-14.78, 16.67), complex(0.75, 0.57)], [complex(-3.70, -21.77), complex(-13.36, 17.55), complex(0.80, 0.47)],
      [complex(-3.62, -21.27), complex(-11.81, 17.93), complex(0.82, 0.44)], [complex(-3.49, -20.66), complex(-10.46, 18.05), complex(0.84, 0.39)],
      [complex(-3.50, -20.29), complex(-9.05, 18.43), complex(0.86, 0.32)], [complex(-3.53, -19.76), complex(-7.65, 18.54), complex(0.88, 0.26)],
      [complex(-3.82, -19.24), complex(-6.01, 18.69), complex(0.89, 0.18)], [complex(-4.15, -19.10), complex(-4.43, 19.05), complex(0.89, 0.10)],
      [complex(-4.18, -18.58), complex(-2.96, 18.78), complex(0.89, 0.04)], [complex(-4.14, -18.03), complex(-1.70, 18.42), complex(0.87, -0.01)],
      [complex(-4.35, -17.39), complex(-0.16, 17.97), complex(0.83, -0.06)], [complex(-4.55, -17.21), complex(1.21, 17.76), complex(0.83, -0.13)],
      [complex(-4.57, -16.62), complex(2.25, 17.15), complex(0.81, -0.16)], [complex(-4.58, -16.05), complex(3.32, 16.44), complex(0.80, -0.21)],
      [complex(-4.86, -16.37), complex(4.79, 16.26), complex(0.77, -0.27)], [complex(-4.60, -15.65), complex(5.57, 15.53), complex(0.74, -0.33)],
      [complex(-4.57, -15.34), complex(6.26, 14.88), complex(0.73, -0.37)], [complex(-3.96, -15.53), complex(6.76, 15.28), complex(0.74, -0.42)],
      [complex(-3.69, -15.41), complex(7.81, 14.72), complex(0.73, -0.46)], [complex(-2.91, -14.52), complex(7.41, 13.68), complex(0.75, -0.46)],
      [complex(-2.66, -12.36), complex(5.84, 12.69), complex(0.76, -0.45)], [complex(-2.51, -12.71), complex(6.62, 12.59), complex(0.73, -0.49)],
      [complex(-2.75, -12.00), complex(7.15, 11.56), complex(0.70, -0.52)], [complex(-2.56, -11.54), complex(7.33, 10.71), complex(0.68, -0.53)],
      [complex(-3.10, -10.94), complex(7.95, 9.41), complex(0.61, -0.54)], [complex(-2.73, -10.60), complex(7.67, 9.08), complex(0.61, -0.57)],
      [complex(-2.96, -10.01), complex(7.99, 7.86), complex(0.58, -0.59)], [complex(-2.22, -9.60), complex(7.68, 7.35), complex(0.57, -0.61)],
      [complex(-2.08, -9.80), complex(8.04, 6.72), complex(0.54, -0.66)], [complex(-2.21, -8.82), complex(7.29, 5.50), complex(0.51, -0.64)],
      [complex(-2.60, -8.83), complex(8.12, 4.66), complex(0.43, -0.66)], [complex(-2.75, -9.06), complex(8.59, 3.79), complex(0.34, -0.71)],
      [complex(-2.77, -8.79), complex(8.49, 3.12), complex(0.30, -0.73)], [complex(-2.18, -8.96), complex(8.49, 3.11), complex(0.26, -0.77)],
      [complex(-2.53, -8.63), complex(8.16, 1.86), complex(0.23, -0.76)], [complex(-1.80, -8.79), complex(8.13, 1.66), complex(0.20, -0.80)],
      [complex(-2.32, -8.12), complex(7.54, 0.56), complex(0.16, -0.74)], [complex(-2.21, -8.04), complex(7.48, -0.19), complex(0.09, -0.75)],
      [complex(-2.90, -8.10), complex(7.16, -1.13), complex(0.02, -0.77)], [complex(-2.67, -8.44), complex(7.32, -1.49), complex(-0.05, -0.78)],
      [complex(-3.26, -7.92), complex(6.73, -2.51), complex(-0.14, -0.73)], [complex(-3.13, -7.85), complex(6.23, -2.75), complex(-0.12, -0.74)],
      [complex(-3.11, -6.58), complex(4.86, -3.03), complex(-0.10, -0.69)], [complex(-2.11, -5.47), complex(4.03, -2.18), complex(-0.06, -0.78)],
      [complex(-1.98, -5.12), complex(3.60, -2.08), complex(-0.10, -0.82)], [complex(-0.97, -5.45), complex(4.32, -1.63), complex(-0.21, -0.85)],
      [complex(-1.02, -4.54), complex(3.50, -1.62), complex(-0.26, -0.84)], [complex(-0.81, -3.71), complex(2.74, -1.30), complex(-0.24, -0.77)],
      [complex(-1.17, -3.22), complex(2.29, -1.49), complex(-0.32, -0.67)], [complex(-0.92, -3.49), complex(2.33, -1.78), complex(-0.32, -0.63)],
      [complex(-0.95, -1.81), complex(1.33, -0.62), complex(-0.33, -0.64)], [complex(0.27, 0.51), complex(0.90, 1.80), complex(-0.23, -0.74)],
      [complex(2.17, 2.33), complex(1.73, 4.05), complex(-0.23, -0.77)], [complex(3.86, 1.32), complex(3.72, 3.78), complex(-0.36, -0.69)],
      [complex(7.39, 0.86), complex(7.16, 4.60), complex(-0.45, -0.74)], [complex(9.18, -0.80), complex(9.73, 3.47), complex(-0.56, -0.77)],
      [complex(7.12, -3.85), complex(8.87, -0.35), complex(-0.44, -0.59)], [complex(4.76, -4.82), complex(7.14, -2.02), complex(-0.55, -0.56)],
      [complex(7.00, -1.61), complex(8.43, 1.01), complex(-0.48, -0.62)], [complex(0.46, -0.26), complex(1.94, 0.49), complex(-0.41, -0.48)],
      [complex(-1.65, -1.43), complex(0.00, -1.11), complex(-0.39, -0.45)], [complex(-0.91, -7.67), complex(1.37, -7.28), complex(-0.50, -0.16)],
      [complex(-2.19, -3.73), complex(-0.34, -3.37), complex(-0.65, -0.25)]]


# Bucle principal
while True:
    if uart.any() >= 8:
        buffer = uart.read(4)
        [re] = struct.unpack('f', buffer)
        if DEBUG:
            print(re, end='')
            print(' ', end='')
        buffer = uart.read(4)
        [im] = struct.unpack('f', buffer)
        if DEBUG:
            print(im)
        
        reb = struct.pack('f', re)
        imb = struct.pack('f', im)
        
        if reb + imb == b'\r\n\r\n\r\n\r\n':
            if DEBUG:
                print('\nFin barrido\n')
                print(measured)
                print()
            i = 0
            medido = medido + 1
            
            if medido == 2:    # Me quedo con el segundo barrido
                medido = 0
                print("S11")
                print(measured)
                print("Permitividad")
                print(permitividad_re_im)
                
                # Guardar los datos en un .txt
                f_s11 = open('s11_' + num_fichero + '.txt', 'w')
                for k in range(len(measured)):
                    f_s11.write(str(measured[k][0]) + '\t' + str(measured[k][1]) + '\n')
                f_s11.close()
                
                f_perm = open('perm_' + num_fichero + '.txt', 'w')
                for k in range(len(permitividad_re_im)):
                    f_perm.write(str(permitividad_re_im[k][0]) + '\t' + str(permitividad_re_im[k][1]) + '\n')
                f_perm.close()
                
                num_fichero = str(int(num_fichero) + 1)
                
                for i in range(len(mensajeOFF)):
                    if DEBUG:
                        print(mensajeOFF[i], end='')
                    uart.write(str(mensajeOFF[i]))
                    if DEBUG:
                        time.sleep(0.1)
        
        else:
            # Relleno de matriz con los valores de S11
            measured[i][0] = re
            measured[i][1] = im
            
            # Relleno de matriz con los valores de la permitividad
            S11 = complex(re, im)
            permitividad_compleja[i] = (c[i][0]*S11 - c[i][1])/(c[i][2] - S11)
            r, phi = cmath.polar(permitividad_compleja[i])
            real, imag = r*math.cos(phi), r*math.sin(phi)
            permitividad_re_im[i][0] = real
            permitividad_re_im[i][1] = imag
            
            if DEBUG:
                print(measured)
            
            i = i + 1
            
            
    # SDI-12
    direccion_sdi12.value(0)
    if uart_sdi12.any():
        ch_sdi12 = uart_sdi12.read(1)
        buffer_sdi12 += ch_sdi12
        if DEBUG:
            print('Buffer: ', end='')
            print(buffer_sdi12)
        if ch_sdi12 == b'!':
            direccion_sdi12.value(1)
            if b'?!' in buffer_sdi12: # ?!
                uart_sdi12.write(ID + '\r\n')
                time.sleep(0.04)
            
            elif bytes(ID, 'utf-8') + b'I!' in buffer_sdi12: # aI!
                uart_sdi12.write(ID + '14jgarrido000001001VNA\r\n')
                time.sleep(0.26)
            
            elif bytes(ID, 'utf-8') + b'A' in buffer_sdi12: # aAb!
                index = buffer_sdi12.index(b'A')
                ID = chr(buffer_sdi12[index + 1])
                uart_sdi12.write(ID + '\r\n')
                time.sleep(0.04)
            
            elif bytes(ID, 'utf-8') + b'M!' in buffer_sdi12: # aM!
                uart_sdi12.write(ID + '0011\r\n')
                time.sleep(0.08)
            
            elif bytes(ID, 'utf-8') + b'D0!' in buffer_sdi12: # aD0! -> S11
                mensaje = ID
                for i in range(n_puntos):
                    if measured[i][0] >= 0:
                        mensaje += '+'
                        mensaje += str(measured[i][0])
                    else:
                        mensaje += str(measured[i][0])
                    if measured[i][1] >= 0:
                        mensaje += '+'
                        mensaje += str(measured[i][1])
                    else:
                        mensaje += str(measured[i][1])
                mensaje += '\r\n'
                uart_sdi12.write(mensaje)
                time.sleep(0.30)
            
            elif bytes(ID, 'utf-8') + b'D1!' in buffer_sdi12: # aD1! -> Permitividad
                mensaje = ID
                for i in range(n_puntos):
                    if permitividad_re_im[i][0] >= 0:
                        mensaje += '+'
                        mensaje += str(permitividad_re_im[i][0])
                    else:
                        mensaje += str(permitividad_re_im[i][0])
                    if permitividad_re_im[i][1] >= 0:
                        mensaje += '+'
                        mensaje += str(permitividad_re_im[i][1])
                    else:
                        mensaje += str(permitividad_re_im[i][1])
                mensaje += '\r\n'
                uart_sdi12.write(mensaje)
                time.sleep(0.30)
                
            elif bytes(ID, 'utf-8') + b'!' in buffer_sdi12: # a!
                uart_sdi12.write(ID + '\r\n')
                time.sleep(0.04)
            
            buffer_sdi12 = b''