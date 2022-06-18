# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 17:15:26 2022

@author: JosephVermeil

https://docs.opencv.org/3.4/db/d5b/tutorial_py_mouse_handling.html
"""

# %%

import cv2 as cv

events = [i for i in dir(cv) if 'EVENT' in i]
print( events )

# %%

import numpy as np
import cv2 as cv

# mouse callback function
def draw_circle(event,x,y,flags,param):
    if event == cv.EVENT_LBUTTONDBLCLK:
        cv.circle(img,(x,y),100,(255,0,0),-1)
        
        
# Create a black image, a window and bind the function to window
img = np.zeros((512,512,3), np.uint8)

cv.namedWindow('image')

cv.setMouseCallback('image',draw_circle)

while(1):
    cv.imshow('image',img)
    if cv.waitKey(20) & 0xFF == 27: # escape
        break
    
cv.destroyAllWindows()

# %%

import numpy as np
import cv2 as cv
drawing = False # true if mouse is pressed
mode = True # if True, draw rectangle. Press 'm' to toggle to curve
ix,iy = -1,-1
# mouse callback function
def draw_circle(event,x,y,flags,param):
    global ix,iy,drawing,mode,img,img_copy
    if event == cv.EVENT_LBUTTONDOWN:
        drawing = True
        ix,iy = x,y
        img_copy = np.copy(img)
    elif event == cv.EVENT_MOUSEMOVE:
        if drawing == True:
            if mode == True:
                img = np.copy(img_copy)
                cv.rectangle(img,(ix,iy),(x,y),(0,255,0),-1)
            else:
                cv.circle(img,(x,y),5,(0,0,255),-1)
    elif event == cv.EVENT_LBUTTONUP:
        drawing = False
        if mode == True:
            cv.rectangle(img,(ix,iy),(x,y),(0,255,0),-1)
        else:
            cv.circle(img,(x,y),5,(0,0,255),-1)

img = np.zeros((252,252,3), np.uint8)
img_copy = np.copy(img)
cv.namedWindow('image')
cv.setMouseCallback('image',draw_circle)
i = 0
while(1):
    cv.imshow('image',img)
    k = cv.waitKey(1) & 0xFF
    if k == ord('m'):
        mode = not mode
    elif k == 27:
        break
cv.destroyAllWindows()

# %%

img = np.zeros((3,3,3), np.uint8)
img_copy1 = img[:][:][:]
img_copy2 = np.copy(img)
img[1][1] = (50,50,50)
print(img_copy1)
print(img_copy2)

