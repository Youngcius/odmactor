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

    def get_data_with_time(self, interval=1e-6, num=100):
        self.cache.clear()
        for _ in range(num):
            time.sleep(interval)
            self.cache.append(self.magnitude)
        return self.cache


# ============================
import asyncio
import numpy as np
from ipywidgets import Button
import plotly.graph_objs as go

N = 1000
lockin = LockInAmplifier()

# create a figure widget and a plot
fig_trace = go.FigureWidget()
# fig_trace.add_scatter(x=trace.getIndex(), y=trace.getData()[0])
fig_trace.add_scatter(x=range(N), y=lockin.get_data_with_time(num=N))


async def update_trace():
    """Update the plot every 0.1 s"""
    while True:
        # fig_trace.data[0].y = trace.getData()[0]
        fig_trace.data[0].y = lockin.get_data_with_time(num=N)

        await asyncio.sleep(0.1)


# If this cell is re-excecuted and there was a previous task, stop it first to avoid a dead daemon
try:
    task_trace.cancel()
except:
    pass

loop = asyncio.get_event_loop()
task_trace = loop.create_task(update_trace())

# create a stop button
button_trace_stop = Button(description='stop')
button_trace_stop.on_click(lambda a: task_trace.cancel())

display(fig_trace, button_trace_stop)

#
# cache = deque()
#
# fig, ax = plt.subplots()
# # The time vector is fixed. No need to read it on every iteration.
# x = np.array(hist.getIndex())
# line, = ax.plot(x, x * 0)
# ax.set_xlabel('Time (ps)')
# ax.set_ylabel('Counts')
# ax.set_title('Correlation histogram via Pyro-RPC')
# while hist.isRunning():
#     y = hist.getData()
#     line.set_ydata(y)
#     ax.set_ylim(np.min(y), np.max(y))
#     plt.pause(0.1)
