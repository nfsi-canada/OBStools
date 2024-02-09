from ObsQA.imports import *
def getLegendTemplate(colors,labels,title_str='Legend',opacity = 0.7,legx='50px',legy='50px'):
    template = """
    {% macro html(this, kwargs) %}

    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>jQuery UI Draggable - Default functionality</title>
    <link rel="stylesheet" href="//code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">

    <script src="https://code.jquery.com/jquery-1.12.4.js"></script>
    <script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>
    
    <script>
    $( function() {
        $( "#maplegend" ).draggable({
                        start: function (event, ui) {
                            $(this).css({
                                right: "auto",
                                top: "auto",
                                bottom: "auto"
                            });
                        }
                    });
    });

    </script>
    </head>
    <body>

    
    <div id='maplegend' class='maplegend' 
        style='position: absolute; bottom: ##LEGY##; left: ##LEGX##; z-index:9999; border:2px solid grey; background-color:rgba(255, 255, 255, 0.8);
        border-radius:6px; padding: 10px; font-size:14px;'>
        
    <div class='legend-title'>##TITLE##</div>
    <div class='legend-scale'>
    <ul class='legend-labels'>
        ##ITEM##
        ##BOTTOM##

    </ul>
    </div>
    </div>
    
    </body>
    </html>

    <style type='text/css'>
    .maplegend .legend-title {
        text-align: left;
        margin-bottom: 5px;
        font-weight: bold;
        font-size: 90%;
        }
    .maplegend .legend-scale ul {
        margin: 0;
        margin-bottom: 5px;
        padding: 0;
        float: left;
        list-style: none;
        }
    .maplegend .legend-scale ul li {
        font-size: 80%;
        list-style: none;
        margin-left: 0;
        line-height: 18px;
        margin-bottom: 2px;
        }
    .maplegend ul.legend-labels li span {
        display: block;
        float: left;
        height: 16px;
        width: 30px;
        margin-right: 5px;
        margin-left: 0;
        border: 1px solid #999;
        }
    .maplegend .legend-source {
        font-size: 80%;
        color: #777;
        clear: both;
        }
    .maplegend a {
        color: #777;
        }
    </style>
    {% endmacro %}"""

    template = template.replace('##TITLE##',title_str)
    template = template.replace('##LEGX##',legx)
    template = template.replace('##LEGY##',legy)
    for c,l in zip(colors,labels):
        item = "<li><span style='background:" + str(c) + ";opacity:" + str(opacity) + ";'></span>" + str(l) + "</li>"
        template = template.replace('##ITEM##',item)
        template = template.replace('##BOTTOM##',"##ITEM##\n    ##BOTTOM##")
    template = template.replace('    ##ITEM##','')
    template = template.replace('    ##BOTTOM##','')
    
    macro = MacroElement()
    macro._template = Template(template)
    return macro