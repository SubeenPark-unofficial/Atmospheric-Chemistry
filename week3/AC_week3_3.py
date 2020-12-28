import datetime as dt

#1
print ("="*15 + " Ex1 " + "="*15)
rightnow = dt.datetime.now()
print (rightnow.strftime("%Y%m%d%H%M%S"))


#2
print ("="*15 + " Ex2 " + "="*15)
time = []
start = dt.datetime(2016,5,1)
numdays = 61

for x in range (0, numdays):
    add = start + dt.timedelta(days = x)
    time.append(add.strftime('%Y-%m-%d'))

print(time)
