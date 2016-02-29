def sqrtit(a):
    x_old = -1
    x = a
    while x_old != x:
        x_old = x
        x = 0.5 * (x + a / x)
    return x


def example():
    print('sqrtit(2) = {}'.format(sqrtit(2)))


if __name__ == '__main__':
    example()
