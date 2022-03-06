import threading
import time
from functools import partial


def acquire(i):
    time.sleep(2)
    print('第{}次loop中采集数据，当前线程：{}'.format(i, threading.currentThread().name))
    time.sleep(2)


if __name__ == '__main__':
    n = 9
    for i in range(1, n + 1):
        print('第{}次loop'.format(i))
        t = threading.Thread(target=partial(acquire, i), name='thread {}'.format(i))
        t.start()
        # time.sleep(1)
        # t.join()
        print('当前线程:{}'.format(threading.currentThread().name))
    print('现在是线程{}'.format(threading.currentThread().name))
