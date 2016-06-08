from Tkinter import *
from tkFont import *






class AutoListbox(Listbox):
    def autowidth(self, maxwidth, list=None):
        f = Font(font=self.cget("font"))
        pixels = 0
        if not list:
            for item in self.get(0, "end"):
                pixels = max(pixels, f.measure(item))
        else:
            for item in list:
                pixels = max(pixels, f.measure(item))

        # bump listbox size until all entries fit
        pixels = pixels + 10
        width = int(self.cget("width"))
        for w in range(0, maxwidth + 1, 5):
            if self.winfo_reqwidth() >= pixels:
                break
            self.config(width=width + w)
