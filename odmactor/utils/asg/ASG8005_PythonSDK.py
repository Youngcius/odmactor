import os
from ctypes import *
import platform

STATUS_CALLBACK = CFUNCTYPE(None, c_int, c_char_p)
STATUS_CALLBACK_COUNT = CFUNCTYPE(None, c_int, c_int, POINTER(c_uint32))


class ASG8005:
    _instance = None

    m_CountCount = 0
    py_callback = {}

    def __new__(cls, *args, **kw):
        if cls._instance is None:
            cls._instance = object.__new__(cls, *args, **kw)
        return cls._instance

    def __init__(self):
        wd = os.path.abspath(os.path.dirname(__file__))
        arch = platform.architecture()[0]

        dll_path = ""

        if arch == '64bit':
            dll_path = os.path.join(wd, 'ASGDLL_x64.dll')
            print(" USE ASGDLL_x64.dll ")
        else:
            dll_path = os.path.join(wd, 'ASGDLL_x86.dll')
            print(" USE ASGDLL_x86.dll ")

        if os.path.isfile(dll_path):
            self.__dll = CDLL(dll_path)
        else:
            raise Exception("can not found dll")

        # dev
        self.__dll.open.restype = c_int
        self.__dll.close_usb.restype = c_int
        self.__dll.monitorDeviceStatus.restype = c_int

        # pub
        self.__dll.setCallbackFunc.restype = c_int
        self.__dll.setCallbackFunc_int.restype = c_int
        self.__dll.getDllInfomation.restype = c_char_p
        self.__dll.start_download.restype = c_int
        self.__dll.stop_download.restype = c_int

        # asg
        self.__dll.pulse_download.argtypes = [POINTER(POINTER(c_double)), POINTER(c_int)]
        self.__dll.pulse_download.restype = c_int
        self.__dll.trigger_download.restype = c_int

        # count
        self.__dll.set_counter_repeat.restype = c_int
        self.__dll.set_counter_repeat.argtypes = [c_int]
        self.__dll.isCountContinu.restype = c_int
        self.__dll.isCountContinu.argtypes = [c_int]
        self.__dll.countTimeStep.restype = c_int
        self.__dll.countTimeStep.argtypes = [c_int]
        self.__dll.countConfig.argtypes = [c_int, c_int]
        self.__dll.counter_download.restype = c_int
        self.__dll.counter_download.argtypes = [POINTER(c_int), c_int]

    # dev
    def connect(self):
        return self.__dll.open()

    def close_device(self):
        return self.__dll.close_usb()

    def get_monitor_status(self):
        return self.__dll.monitorDeviceStatus()

    # pub
    def set_callback(self, func):
        if type(func) == STATUS_CALLBACK:
            return self.__dll.setCallbackFunc(func)
        else:
            return False

    def set_callback_count(self, func):
        if type(func) == STATUS_CALLBACK_COUNT:
            return self.__dll.setCallbackFunc_int(func)
        else:
            return False

    def get_device_info(self):
        return str(self.__dll.getDllInfomation())

    def start(self, count=1):
        return self.__dll.start_download(count)

    def stop(self):
        return self.__dll.stop_download()

    # asg
    def checkdata(self, asg_data, length):
        channelLen = [0, 0, 0, 0, 0, 0, 0, 0]
        for i in length:
            if i % 2 != 0 or i < 2:
                return bool(False)
        for i in range(len(asg_data)):
            if len(asg_data[i]) != length[i]:
                return bool(False)
            if len(asg_data[i]) == 2:
                if ((asg_data[i][0] < 7.5) and (asg_data[i][0] != 0)) or (asg_data[i][0] > 26000000000.0) \
                        or ((asg_data[i][1] < 10) and (asg_data[i][1] != 0)) or (asg_data[i][1] > 26000000000.0):
                    return bool(False)
                continue
            for j in range(0, len(asg_data[i]) - 1, 2):

                aint = int(asg_data[i][j] * 1000000)
                bint = int(asg_data[i][j + 1] * 1000000)

                afloat = int(asg_data[i][j] * 100) * 10000
                bfloat = int(asg_data[i][j + 1] * 100) * 10000

                if (aint != afloat or bint != bfloat) or (aint % 50000 != 0 or bint % 50000 != 0):
                    return bool(False)

                if j == 0:
                    if ((asg_data[i][0] < 7.5) and (asg_data[i][0] != 0)) or (asg_data[i][0] > 26000000000.0) or (
                            asg_data[i][1] < 10) or (asg_data[i][1] > 26000000000.0):
                        return bool(False)
                elif j == len(asg_data[i]) - 2:
                    if (asg_data[i][j] < 7.5) or (asg_data[i][j] > 26000000000.0) or (
                            (asg_data[i][j + 1] < 10) and (asg_data[i][j + 1] != 0) or (
                            asg_data[i][j + 1] > 26000000000.0)):
                        return bool(False)
                else:
                    if (asg_data[i][j] < 7.5) or (asg_data[i][j] > 26000000000.0) or (asg_data[i][j + 1] < 10) or (
                            asg_data[i][j + 1] > 26000000000.0):
                        return bool(False)
                channelLen[i] += (asg_data[i][j] + asg_data[i][j + 1])
        for i in range(8):
            if channelLen[i] > 5200000000:
                return bool(False)
        return bool(True)

    def download_ASG_pulse_data(self, asg_data, length):
        if True != self.checkdata(asg_data, length):
            exit(" ASG Data  error !")
        c_length = (c_int * 8)(*tuple(length))
        max = 0
        for i in range(8):
            if max < length[i]:
                max = length[i]
        c_asg_data = (c_double * max * 8)(*(tuple(i) for i in asg_data))
        c_asg_data = (POINTER(c_double) * len(c_asg_data))(*c_asg_data)
        return self.__dll.pulse_download(c_asg_data, c_length)

    def ASG_trigger_download(self):
        return self.__dll.trigger_download()

    # count
    def ASG_set_counter_repeat(self, repeat):
        repeat = repeat * 2
        return self.__dll.set_counter_repeat(c_int(repeat))

    def ASG_isCountContinu(self, isContinu):
        return self.__dll.isCountContinu(c_int(isContinu))

    def ASG_countTimeStep(self, timeStep):
        return self.__dll.countTimeStep(c_int(timeStep))

    def ASG_countConfig(self, isCountEnable, asgConfig=0xff):
        return self.__dll.countConfig(c_int(asgConfig), c_int(isCountEnable))

    def checkCountData(self, countData, length):
        countLength = 0
        for i in range(0, length, 2):
            if (countData[i] < 20) or (countData[i] % 5 != 0) or (countData[i] != int(countData[i])):
                return bool(False)
            if (countData[i + 1] < 5) or (countData[i + 1] % 5 != 0) or (countData[i + 1] != int(countData[i + 1])):
                return bool(False)
            countLength += (countData[i] + countData[i + 1])
        if countLength < 1500:
            return bool(False)
        return bool(True)

    def ASG_counter_download(self, count_data, length):
        if True != self.checkCountData(count_data, length):
            exit(" Count Data  error !")
        m_CountCount = 1
        count_data = (c_int * len(count_data))(*tuple(count_data))
        return self.__dll.counter_download(count_data, length)
