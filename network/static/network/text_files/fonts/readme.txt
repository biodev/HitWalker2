
This zip file contains a set of regular and bold Ubuntu fonts, together with some CSS code for embedding them into your web page.

The handling of fonts (and other assets) is different across KeyLines Flash & HTML5 components.

In the case of Flash, the fonts and other assets are already embedded into the .swf file.

In HTML5, fonts and other assets are loaded like any other resource in the web page.  KeyLines looks for the 'Ubuntu' font.  This font must be already loaded into your webpage by the time KeyLines creates a chart.  In practical terms this means waiting for the $(window).load function before using KeyLines: see the example code for how to do this.

If the text does not show immediately, the chances are that the font has not loaded into the browser page by the time the chart is created: using $(window).load will guarantee KeyLines has the resources it needs.

To use the contents of this zip file:

1) Unpack all the files to a directory served by your web server
2) Take the contents of the 'stylesheet.css' file and amend the paths so they point to the correct directory - then put this into your own .css file for your web pages
3) Test that the files are accessible

