from Tkinter import *
import tkMessageBox
def antwort():
    lab=Label(root,text="Hier nicht! ")  
    lab.pack()
    w = Canvas(root, width=200, height=100,bg='red')
    w.create_rectangle(50, 25, 150, 75, fill="blue")
    w.pack()
    #tkMessageBox.showinfo('Hier nicht!','Hier auch nicht!')
root=Tk()
but=Button(root,text="Wo ist Tommy?",command=antwort)
but.pack()
root.mainloop()