%.py: %.ui
	pyuic $< -o $@

UI_GENERATED = dartnavigator.py displaysettings.py

all: $(UI_GENERATED)

clean:
	rm -f $(UI_GENERATED)
