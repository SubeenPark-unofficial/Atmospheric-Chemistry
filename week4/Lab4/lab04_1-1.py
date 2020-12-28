# 실습 1-1: 피보나치 수열 구하기

# 양의 정수 n 입력받기
n = int(input("how many terms? "))

# 변수 초기화
x = 0
y = 1

# TODO. while loop을 사용해 피보나치 수열 출력하기
while n > 0:
    print(x, end = " ")
    temp = x
    x = y
    y += temp
    n -= 1
