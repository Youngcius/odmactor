import time
import TimeTagger as tt

tagger = tt.createTimeTagger()
tagger.setTestSignal(1,True) # ~1 MHz
tagger.setTestSignal(2,True)

# cnt = tt.CountBetweenMarkers(tagger, 2, 1, -1, n_values=10000)
cnt = tt.Counter(tagger, channels=[1,2], binwidth=int(1e9), n_values=int(1e3))
cnt.startFor(int(0.5e12)) # 0.5 s
# time.sleep(0.5)
time.sleep(3)

data = cnt.getData()
print(data.shape)
print(data)
print(data.sum())

print(cnt.isRunning())
tt.freeTimeTagger(tagger)
