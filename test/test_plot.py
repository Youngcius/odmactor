# import matplotlib
# gui_env = ['TKAgg','GTKAgg','Qt4Agg','WXAgg']
# for gui in gui_env:
#     try:
#         print("testing", gui)
#         matplotlib.use(gui, force=True)
#         from matplotlib import pyplot as plt
#         break
#     except:
#         continue
# print("Using:",matplotlib.get_backend())


import matplotlib
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
print(plt.get_backend())

plt.plot([1, 2, 3])
plt.show()



# gui_env = [i for i in matplotlib.rcsetup.interactive_bk]
# non_gui_backends = matplotlib.rcsetup.non_interactive_bk
# print("Non Gui backends are:", non_gui_backends)
# print("Gui backends I will test for", gui_env)
# for gui in gui_env:
#     print("testing", gui)
#     try:
#         matplotlib.use(gui, warn=False, force=True)
#         from matplotlib import pyplot as plt
#
#         print("    ", gui, "Is Available")
#         plt.plot([1.5, 2.0, 2.5])
#         fig = plt.gcf()
#         fig.suptitle(gui)
#         plt.show()
#         print("Using ..... ", matplotlib.get_backend())
#     except:
#         print("    ", gui, "Not found")
