class ProgressHook(object):
    """Represents a sub-progress range.  When called with a number
    between 0.0 and 1.0, calls the next hook with the scaled range.
    Sub-ranges can be created with subProgress."""
    
    def __init__(self, hook, range = (0.0, 1.0)):
        self._start = range[0]
        self._end = range[1]
        self._hook = hook

    def start(self):
        return self._start

    def end(self):
        return self._end

    def __call__(self, progress):
        self._hook(self._start + progress * self._length())

    def _length(self):
        return self._end - self._start

    def rangeTicker(self, count):
        return RangeTicker(self, count)

    def subProgress(self, startPercent, endPercent):
        return ProgressHook(
            self._hook,
            (self.start() + startPercent * self._length(),
             self.start() + endPercent * self._length()))

class RangeTicker(object):
    def __init__(self, progressHook, count):
        self._progressHook = progressHook
        self._count = float(count)
        self._index = 0

    def __call__(self):
        self._index += 1
        self._progressHook(self._index / self._count)

# --------------------------------------------------------------------

import sys, time

class StatusMessage(object):
    def __init__(self, message, stream = None):
        self._message = message
        self._promille = 0
        self._startClock = time.clock()
        if stream is None:
            stream = sys.stdout
        self._stream = stream
        self._stream.write("\r%s... " % (self._message, ))
        self._stream.flush()

    def __del__(self):
        self.finish()

    def __call__(self, percent):
        promille = int(percent * 1000)
        if promille != self._promille:
            self._stream.write("\r%s... (%.1f%%)" % (
                self._message, percent * 100))
            self._stream.flush()
            self._promille = promille

    def totalTime(self):
        return self._totalTime

    def finish(self):
        if self._promille is not None:
            self._totalTime = time.clock() - self._startClock
            self._stream.write("\r%s... done. (%ss.)\n" % (
                self._message, self._totalTime))
            self._promille = None

if __name__ == "__main__":
    sm = StatusMessage("please wait")
    p1 = ProgressHook(sm)
    p1(0.1)
    p1(0.2)
    p1(0.45)
    p2 = p1.subProgress(0.5, 1.0)
    p3 = p2.subProgress(0, 0.5)
    p3(0.0)
    p3(0.5)
    p3(0.99)
    p2(0.8)
    p2(1.0)
    sm.finish()
