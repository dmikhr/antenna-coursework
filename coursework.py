# -*- encoding: utf-8 -*-

# Исследование антенны бегущей волны с замедленной переменной фазовой скоростью

# импорт библиотек, необходимых для проведения вычислений
import os
import numpy as np
from scipy.integrate import quad
from scipy import real, imag
from numpy import cos, sin, pi, exp, abs
from pylab import *

#os.chdir(os.path.dirname(__file__))

# функция для интегрирования комплексных функций
def compl_integr(f, a, b):
 # функция quad возвращает резульат интегрирования в виде вектора, содерж. 2 элемента: 
 # 0-ой - результат интегрирования, 1-ый - точность вычисления
 RE = quad( lambda x: real( f(x) ), a, b )[0] # интегрирование действительной части
 IM = quad( lambda x: imag( f(x) ), a, b )[0] # интегрирование мнимой части
 return RE + 1j*IM # реузльтат в виде комплексного числа

# функция для построения графиков
def grafik(X, Y, fname, variant, dlina_antenni, limit=-1, sizex=10, sizey=5):
 figure(figsize=(sizex, sizey)) # размеры графика
 ax = subplot(111) # инициализация
 # задание границы графика
 if limit > 0:
  ax.set_xlim(right=limit)
 ax.plot(X, Y)
 # добавление сетки на график
 ax.grid()
 # путь, по которому будет сохранён файл
 savepath = 'plots/{0}/{1}_{2}.png'.format(str(dlina_antenni), variant, fname)
 savefig( savepath ,format='PNG') # сохранение графика в виде изображения

# нормировка
def normir(data):
 M = np.amax(data) # поиск максимального значения
 # каждый элемент массива делится на максимальное число в массиве
 for id in xrange(len(data)): data[id] =  data[id] / M
 return data # возврат результата нормировки

# длина волны 0.03 м
lambd0 = 0.03

# расчёт K0
K0 = (2 * pi) / float(lambd0)

# длина антенны
L = lambda i: i*lambd0

# кси оптимальное
ksi_opt = lambda i: 1 + lambd0 / float(2*L(i))

# постоянная распротранения Г оптимальная
G_opt = lambda i: K0 * ksi_opt(i)

# расчёт delta
delta = lambda i: G_opt(i) - K0

# словарь функций постоянной распространения Г (4 варианта)
Gammas = {'a': lambda z, n, i: G_opt(i) + n * (delta(i) / 10)*(1 - z/L(i)),\
    'b': lambda z, n, i: G_opt(i) + n * (delta(i) / 10)*(1 - 2 * z/L(i)),\
    'c': lambda z, n, i: G_opt(i) - (z/L(i)) * (n * (delta(i) / 10)),\
    'd': lambda z, n, i: G_opt(i) + (z/L(i)) * (n * (delta(i) / 10))}

# расчёт диаграммы направленности
def f( theta, n, i , Gfunc):
 # подынтегральное выражение
 fn = lambda z: exp( 1j * K0 * z * cos( theta )  - 1j * Gfunc(z, n, i) * z ) 
 compl_res =  compl_integr(fn, 0, L(i)) # расчёт интеграла
 return abs(compl_res) # модуль результата интегрирования 

# расчёт коэффициента направленного действия
def D(n, i , Gfunc):
 # подынтегральное выражение
 func = lambda theta: ((f(theta, n, i , Gfunc))**2)*sin(theta)
 # интегрирование
 compl_res =  compl_integr(func, 0 , pi)
 # результат, как D=2/(результат интегрирования)
 return (2/abs(compl_res))

# задание массива n=0,2,4,6,8,10
nrange = xrange(0, 11, 2)
# задание расчётного интервала с шагом для f(theta) 
xitems = np.arange(0, 3.2, 0.01)
xitems_len = len(xitems)
# инициализация массва для результатов расчёта диаграммы направленности
yitems = [np.zeros(xitems_len) for i in xrange(len(nrange))]

# задание расчётного интервала с шагом для КНД
dxitems = np.arange(0, 11, 2)
# инициализация массва для результатов расчёта КНД
dyitems = np.zeros(len(dxitems))

# инициализация массивов для расчёта theta0.7 и УБЛ1
y_2theta07 = np.zeros(len(dxitems))
UBL1 = np.zeros(len(dxitems))

# длина антенн: 3, 6 и 10 длин волн
#dlina_antenn = [3, 6, 10]
dlina_antenn = [3]

for i in dlina_antenn: # произвести расчёт для каждой длины антенны
 print 'antenna', i, 'lambda'
 print 'delta', delta(i)
 for variant, Gamma in Gammas.items(): # расчёт для каждой Г
  # нормировка
  print 'variant', variant
  dc = 0 # инициализация счётчика отсчётов для КНД 
  # расчёт для n=0,2,4,6,8,10
  for n in nrange:
   M = f(0, n, i, Gamma) # нормировочный коэффициент
   c = 0 # инициализация счётчика отсчётов для f(theta)
   dyitems[dc] = D(n, i, Gamma) # расчёт КНД
   dc += 1 #  инкремент счётчика для КНД 
   theta_07_found = False # theta0.7 пока не найденна
   UBL1_found = False # УБЛ1 пока не найден
    # расчёт f(theta) для theta в заданном диапазоне
   for theta in xitems:
    # вычисление f(theta)
    yitems[n/2][c] = f(theta, n, i, Gamma)/M # расчёт диаграммы направленности
    # определение значения 2theta0.7
    if theta_07_found == False: # если theta0.7 ещё не найденна
     #если текущее значение f(theta) стало меньше 0.7 
     if yitems[n/2][c] < 0.7: 
      # то в качестве theta0.7 сохраняем предыдущее значение f(theta)
      if n == 0:
       theta07_opt = tmp_theta
      y_2theta07[n/2] = tmp_theta / theta07_opt
      theta_07_found = True # theta0.7 найденна
     else:
      # если f(theta) пока болше 0.7, то сохр. значение в доп. переменную
      tmp_theta = theta 
    # определение УБЛ1
    if UBL1_found == False: # если УБЛ1 ещё не найден
     if c > 1: # чтобы найти экстремум, должно пройти не менее 3-х отсчётов
      # поиск локального экстремума
      if yitems[n/2][c - 2] < yitems[n/2][c - 1] > yitems[n/2][c]:
       # если экстремум найден, сохр его как УБЛ1 
       if n == 0:
        UBL1_opt = yitems[0][c - 1] 
       UBL1[n/2] = yitems[n/2][c - 1] / UBL1_opt 
       UBL1_found = True # УБЛ1 для заданного n найден
      
    c += 1 # инкремент счётчика отсчётов для f(theta)
   print 'sdelano', n
  
  #нормировка для D
  dyitems = normir(dyitems)
  
  # построение графика диаграммы направленности
  figure(figsize=(10,5)) # размер графика
  ax = subplot(111) # инициализация нового графика
  ax.set_xlim(right=pi) # предел по оси x - пи
  ax.grid() # добавление сетки
  for n_id in nrange: # построение графика f(theta) для каждого n
   yitems[n_id/2] = normir( yitems[n_id/2] )
   ax.plot(xitems, yitems[n_id/2]) 
  # сохранение в файл
  savefig( 'plots/'+str(i)+'/'+variant+'_plot.png', format='PNG' )
  
  # построение графика КНД(n)
  grafik(dxitems, dyitems, 'D', variant, i) 
  # построение графика 2theta07(n)
  grafik(dxitems, y_2theta07, '2theta07', variant, i)
  # построение графика для УБЛ1(n)  
  grafik(dxitems, UBL1, 'ubl1', variant, i)
