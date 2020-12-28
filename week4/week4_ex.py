#import
import pandas as pd
import datetime
import numpy as np
import matplotlib.pyplot as plt
#data setting
data=pd.read_excel('week4_ex.xls', encoding = 'euc.kr')
data_num=data.drop([0], axis=0)
data_num2=data_num.rename(columns = {'날짜':'date','오 존': 'O3','이산화질소':'NO2','일산화탄소':'CO','아황산가스':'SO2'})
finedust=data_num2[['PM10','PM2.5']]
others=data_num2[['O3','NO2','CO','SO2']]


plt.subplots(2,1)

#plot
finedust.plot(style='o-',ms=3)
plt.xlabel('Day')
plt.ylabel('Daily average amount of fine dust '+r'$(μg/m^3)$')
plt.title('Daily average amount of fine dust in February 2019')


others.plot.bar(stacked=True)
plt.xlabel('Day')
plt.title('Daily average amount of air pollutants('+r'$O_3,NO_2,CO,SO_2$'+') in February 2019')
plt.ylabel('Daily average amount of air pollutants(ppm)')

plt.show()
