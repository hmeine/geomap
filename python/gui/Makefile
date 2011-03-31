PYUIC = pyuic4
PYRCC = pyrcc4

%_ui.py: %.ui
	$(PYUIC) $< -o $@

%_rc.py: %.qrc
	$(PYRCC) $< -o $@

UI_GENERATED = dartnavigator_ui.py displaysettings_ui.py icons_rc.py

all: $(UI_GENERATED)

docs:
	epydoc --no-sourcecode --no-frames --no-private -v --docformat restructuredtext --output=/home/meine/public_html/geomap/api geomap.GeoMap delaunay maputils tools statistics levelcontours pixelmap crackConvert dartpath

clean:
	rm -f $(UI_GENERATED)

.PHONY: docs
