from Tkinter import *

class ImageFrame(Frame):
    def __init__(self, _tk_, file_path, **options):
        Frame.__init__(self, _tk_, **options)

        #DesignPhoto =PhotoImage(file = file_path)
        DesignPhoto = PhotoImage(file = file_path)
        self.Photo = Label(master=self, image=DesignPhoto)
        self.Photo.image = DesignPhoto

        self.Photo.grid(row=0, column=0, padx=3, pady=3, sticky=W+E+N+S)


