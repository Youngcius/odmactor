import time
import pyvisa
from pymeasure import instruments


class LockInAmplifier(instruments.srs.SR830):
    def __init__(self, **kwargs):
        rm = pyvisa.ResourceManager()
        rs_list = rm.list_resources()
        found = False
        for rs in rs_list:
            if rs.startswith('GPIB0::') and rs.endswith('INSTR'):
                super(LockInAmplifier, self).__init__(rs, **kwargs)
                print('Lock-in Amplifier is found on {}:'.format(rs), self.id)
                found = True
                break
        if not found:
            print('Lock-in Amplifier is not found!')

        self.cache = []

    def get_data_with_time(self, interval=1e-3, num=100):
        self.cache.clear()
        for _ in range(num):
            time.sleep(interval)
            self.cache.append(self.magnitude)
        return self.cache


cache = deque()

fig, ax = plt.subplots()
# The time vector is fixed. No need to read it on every iteration.
x = np.array(hist.getIndex())
line, = ax.plot(x, x * 0)
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Counts')
ax.set_title('Correlation histogram via Pyro-RPC')
while hist.isRunning():
    y = hist.getData()
    line.set_ydata(y)
    ax.set_ylim(np.min(y), np.max(y))
    plt.pause(0.1)
