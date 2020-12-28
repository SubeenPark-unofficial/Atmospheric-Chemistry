# 실습 1-2: 피보나치 수열 구하기

# 양의 정수 n 입력받기
n = int(input("how many terms? "))

# 변수 초기화
x = 0
y = 1

# TODO. for loop을 사용해 피보나치 수열 출력하기
for i in range(n):
    print(x, end = " ")
    temp = x
    x = y
    y += temp
