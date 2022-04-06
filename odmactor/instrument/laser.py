class Laser(object):
    def __init__(self):
        pass

    def connect(self):
        """
        Reconnect the instrument
        """
        print('Laser 重新连接')

    def start(self):
        print('Laser 打开')

    def stop(self):
        print('Laser 关闭')

    def close(self):
        print('Laser 释放')

    def set_power(self, power):
        print('设置 Laser 功率为', power)
