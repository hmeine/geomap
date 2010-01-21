%.py: %.ui
	pyuic $< -o $@

UI_GENERATED = dartnavigatorbase.py displaysettings.py

all: $(UI_GENERATED)

docs:
	epydoc --no-sourcecode --no-frames --no-private -v --docformat restructuredtext --output=/home/meine/public_html/geomap/api hourglass.GeoMap delaunay maputils tools statistics levelcontours pixelmap crackConvert dartpath

clean:
	rm -f $(UI_GENERATED)

.PHONY: docs
