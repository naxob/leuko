import Tkinter
import Canvas

import math
 
import sys  
from Tkconstants import *

class Butterfly:

  def __init__(self, master = None):
    self.canvas = Tkinter.Canvas(master, relief = RIDGE, 
    bd = 2, bg = "white", width = 300, height = 300)
    self.canvas.pack()
  
    self.button = Tkinter.Button(master, text = " Run ", command = self.draw)
    self.button.pack(side = BOTTOM, pady = 4)

  def draw(self):
    theta  = 0.0
    while theta < 75.39:
      r = math.exp(math.cos(theta))-2*math.cos(4*theta)+(math.sin(theta/12))**5

      x = r*math.cos(theta)
      y = r*math.sin(theta)
 
      xx = (x*30) + 150
      yy = (y*30) + 150
      if (theta == 0.0):
        Canvas.Line(self.canvas, xx, yy, xx, yy)
      else:
        Canvas.Line(self.canvas, xOld, yOld, xx, yy)
      self.canvas.update_idletasks()
      xOld = xx
      yOld = yy
      theta = theta + 0.02
    """      
    def draw(self):
        theta  = 0.0
        while theta < 75.39:
          r = math.exp(math.cos(theta))-2*math.cos(4*theta)+(math.sin(theta/12))**5
    
          x = r*math.cos(theta)
          y = r*math.sin(theta)
     
          xx = (x*30) + 150
          yy = (y*30) + 150
          if (theta == 0.0):
            Canvas.Line(self.canvas, xx, yy, xx, yy)
          else:
            Canvas.Line(self.canvas, xOld, yOld, xx, yy)
          self.canvas.update_idletasks()
          xOld = xx
          yOld = yy
          theta = theta + 0.02
    """
    def exit(self):
            sys.exit(0)

if __name__ == "__main__":

   root = Tkinter.Tk()
   root.title("Butterfly Curve")
 

   butterfly = Butterfly(root)

   root.mainloop()