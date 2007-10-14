import copy, qt
import mapdisplay, icons

class Workspace(mapdisplay.MapDisplay):
    __base = mapdisplay.MapDisplay

    __slots__ = ("_level0", "_mapRestartAction")
    
    def __init__(self, level0, originalImage):
        self.__base.__init__(self, copy.deepcopy(level0), originalImage)
        self._level0 = level0

        reinitIcon = qt.QPixmap()
        reinitIcon.loadFromData(icons.reinitIconPNGData, "PNG")
        ra = qt.QAction(self, "mapRestartAction")
        ra.setIconSet(qt.QIconSet(reinitIcon))
        ra.setText("Re-start from level 0")
        ra.setMenuText("&Re-start from level 0")
        ra.setToolTip("Re-start with level 0 map")
        ra.addTo(self.Tools)
        self.connect(ra, qt.SIGNAL("activated()"), self.restart)
        self.mapRestartAction = ra

    def _detachMapStats(self):
        for a in self.map.__dict__:
            o = getattr(self.map, a)
            if hasattr(o, "detachHooks"):
                print "detaching hooks of", o
                o.detachHooks()

    def restart(self):
        """Restart with level 0/1."""

        self._detachMapStats()
        self.setMap(copy.deepcopy(self._level0))
        
