import threading
import time


# 新线程执行的代码:
def loop():
    print('thread %s is running...' % threading.current_thread().name)
    n = 0
    while n < 5:
        n = n + 1
        print('thread %s >>> %s' % (threading.current_thread().name, n))
        time.sleep(1)
    print('thread %s ended.' % threading.current_thread().name)


if __name__ == '__main__':
    print('thread %s is running...' % threading.current_thread().name)
    t = threading.Thread(target=loop, name='LoopThread')
    print('start before:', time.time_ns())

    t.start()
    print('join before:', time.time_ns())



    # time.sleep(15)
    # 即主线程任务结束之后，进入阻塞状态，一直等待其他的子线程执行结束之后，主线程在终止
    t.join(timeout=3) # join的作用是设置该线程在主线程接下来执行之前执行
    print('thread %s ended.' % threading.current_thread().name)
