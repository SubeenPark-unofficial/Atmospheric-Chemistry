# 실습 3: 소수 찾기

# 양의 정수 입력 받기
n = int(input("enter integer:"))

# TODO. nested loop을 이용하여 n 이하의 소수 찾기
for i in range(1, n + 1):
    if i > 1:
        for j in range(2, i):
            if i % j == 0:
                break
        else:
            print(i)
