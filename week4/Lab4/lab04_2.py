# 실습 2: 369게임 출력하기

# TODO 1. while loop을 사용해 1과 99 사이의 정수 n 입력 받기
# 출력 코드:
# n = int(input("enter integer in range 1 ~ 99: "))
# print("integer must be in the range 1 ~ 99")
while True:
    n = int(input("enter integer in range 1 ~ 99: "))
    if n > 0 and n < 100:
        break
    else:
        print("integer must be in the range 1 ~ 99")


# TODO 2. while loop 또는 for loop을 이용해 369게임 출력하기
for i in range(1, n + 1):
    a = i // 10
    b = i % 10
    if a in (3, 6, 9):
        if b in (3, 6, 9):
            print("**", end = " ")
        else:
            print("*", end = " ")
    else:
        if b in (3, 6, 9):
            print("*", end = " ")
        else:
            print(i, end = " ")

    if i % 10 == 0:
        print("\n", end = "")
