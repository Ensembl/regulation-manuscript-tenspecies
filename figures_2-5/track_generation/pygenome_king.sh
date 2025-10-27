#human
pyGenomeTracks --height 12 --fontSize 12 --plotWidth 50 --trackLabelFraction 0.1 --tracks tracks_human_king.ini --region  chr11:5,201,829-5,307,000 --outFile exports/tracks_human_king.svg
pyGenomeTracks --height 12 --fontSize 12 --plotWidth 50 --trackLabelFraction 0.1 --tracks tracks_human_king_zoom.ini --region  chr11:5,266,216-5,303,059 --outFile exports/tracks_human_king_zoom.svg

# mouse
pyGenomeTracks --height 12 --fontSize 12 --plotWidth 50 --trackLabelFraction 0.1 --tracks tracks_mouse_king.ini --region chr7:103,438,805-103,538,500 --outFile exports/tracks_mouse_king.svg
