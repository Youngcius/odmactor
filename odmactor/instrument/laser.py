class Laser(object):
    def __init__(self):
        print('Laser initialized')
        pass

    def connect(self):
        """
        Reconnect the instrument
        """
        print('Laser reconnected')

    def start(self):
        print('Laser opened')

    def stop(self):
        print('Laser closed')

    def close(self):
        print('Laser released')

    def set_power(self, power):
        print('Set Laser power:', power)
