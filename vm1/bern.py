from tkinter import *
import math

# окно
win = Tk()
win.title("Grafik")
win.geometry("600x600+100+30")
canvas = Canvas(win, width=600, height=600, bg="LightSalmon")
canvas.pack()


# оси
def axes():
    canvas.create_line(0, 300, 600, 300, fill="DarkViolet", width=1, arrow=LAST)
    canvas.create_line(300, 0, 300, 600, fill="DarkViolet", width=1, arrow=FIRST)


# грифик
def Lem():
    a = 3
    fi = 0
    sc = 20
    while fi <= (2 * math.pi):
        ro = math.sqrt((2 * math.sqrt(a)))
        ro = ro * math.cos(2 * fi)
        x = sc * ro * math.cos(fi)
        y = sc * ro * math.sin(fi)
        canvas.create_oval(300 + x, 300 - y, 300 + x, 300 - y, fill="blue")
        fi = fi + math.pi / 1800


win.title()
axes()
Lem()
win.mainloop()