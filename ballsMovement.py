# import modules
import math
import random
from scipy.optimize import brentq

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib import colors as mcolors

    # -------------------------------------------------------------------------
    # -- CLASS Position
    # -------------------------------------------------------------------------
    
class Position:
    def __init__(self,
                 x = 0,
                 y = 0):
        self.x = x
        self.y = y
    
    def updateX(self,x):
        self.x = x
    
    def updateY(self,y):
        self.y = y
    
    def getX(self):
        return self.x
    
    def getY(self):
        return self.y
    
    def __str__(self):
        return "position: x = {:.3f}, y {:.3f}".format(self.x,self.y)
    
    def __repr__(self):
        return "position: x = {:.3f}, y {:.3f}".format(self.x,self.y)
        
    
###############################################################################
    
    # -------------------------------------------------------------------------
    # -- CLASS SKELETON
    # -------------------------------------------------------------------------

class Skeleton:
    """
    Since we are studying the evolution of balls, the skeletons used are points
    in the space.
    
    A Skeleton 
    """
    
    def __init__(self,
                 position = Position(),
                 vx = 1,
                 vy = 1):
        
        self.position = position
        self.vx = vx
        self.vy = vy
    
    def updatePos(self,pos):
        self.position = pos
    
    
    def updateSpeed(self,vx,vy):
        self.vx = vx
        self.vy = vy
    
    def getX(self):
        return self.position.x
    
    def getY(self):
        return self.position.y
    
    def getPos(self):
        return self.position
        
    def getVX(self):
        return self.vx
    
    def getVY(self):
        return self.vy

    def getSpeed(self):
        return self.vx, self.vy

###############################################################################
    
    # -------------------------------------------------------------------------
    # -- CLASS BALL
    # -------------------------------------------------------------------------

class Ball:
    """
    Class used to represent the balls, it stores a skeleton (used to locate the
    object in the space), the mass of the ball, a color to display it, a 
    function of the euclidien space used to compute the isosurface and the 
    isovalue.
    """
    def __init__(self,
                 skeleton = Skeleton(),
                 mass = 0.10,
                 color = 'r',
                 potential = None,
                 isovalue = 0.05,
                 border_size = 10):
        
        self.skeleton = skeleton
        self.mass = mass
        self.color = color
        self.isovalue = isovalue
        self.border_size = border_size        
        self.size = 0
    
    def potential(self,pos):
        return math.sqrt((pos.x-self.getX())**2 + (pos.y - self.getY())**2)

    def isInside(self,pos):        ## convention : potential = 0 in the position of the skeleton and increasing    
        return self.potential(pos) < self.isovalue
    
    def getborder(self,nb_points = 50):
        if nb_points == 0:
            nb_points = self.border_size
        angles = [(2 * math.pi * k) / nb_points for k in range(nb_points)]
        border = []
        for angle in angles:
            
            r = brentq(lambda r: self.potential(Position(r * math.cos(angle)+self.getX(),
                                                        r * math.sin(angle)+self.getY())) - self.isovalue,0,5)
            border.append(Position(r * math.cos(angle)+self.getX(),
                                                     r * math.sin(angle)+self.getY()))
        self.border = border
    
    def isColliding(self,ball):
        for p in self.border:
            if ball.isInside(p):
                return True
        return False
    
    def isBouncing(self,space): #gauche, droite, haut, bas
        left = self.left < space[0] and self.getVX() < 0
        right = self.right > space[1] and self.getVX() > 0
        up = self.up > space[2] and self.getVY() > 0
        down = self.down < space[3] and self.getVY() < 0
        return left or right, up or down
    
    def changeColor(self,col):
        self.color = col
    
    def updatePos(self,pos):
        self.skeleton.updatePos(pos)
        self.getborder()
#        print("len(self.border) = {}".format(len(self.border)))
        self.up = max(self.border, key=lambda p:p.y).y
        self.down = min(self.border, key=lambda p:p.y).y
        self.left = min(self.border, key=lambda p:p.x).x
        self.right = max(self.border, key=lambda p:p.x).x
        self.size = (abs(self.up - self.getY()) + abs(self.down - self.getY()) + abs(self.left - self.getX()) + abs(self.right - self.getX()))/4
#        print("pos = {}".format(self.getPos()))
#        print("l = {:.3f}, r = {:.3f}, u = {:.3f}, d = {:.3f} \n".format(self.left, self.right, self.up, self.down))
    
    def updateSpeed(self,vx,vy):
        self.skeleton.updateSpeed(vx,vy)
    
    def bounceX(self):
        self.skeleton.vx *= -1
    
    def bounceY(self):
        self.skeleton.vy *= -1
    
    def getX(self):
        return self.skeleton.getX()
    
    def getY(self):
        return self.skeleton.getY()
    
    def getPos(self):
        return self.skeleton.getPos()
        
    def getVX(self):
        return self.skeleton.getVX()
    
    def getVY(self):
        return self.skeleton.getVY()

    def getSpeed(self):
        return self.skeleton.getSpeed()
    
    def getMass(self):
        return self.mass
    
    def getInside(self):
        x_space = np.linspace(self.getX()-self.isovalue, self.getX() + self.isovalue, num = 20)
        y_space = np.linspace(self.getY()-self.isovalue, self.getY() + self.isovalue, num = 20)
        inside = []
        for x in x_space:
            for y in y_space:
                p = Position(x,y)
                if self.isInside(p):
                    inside.append(p)
        return inside

################################## !!!!!!!!!!!!!!!!!!!!!!!!!!#################"
        x_space = np.linspace(self.left, self.right, num = 10, endpoint=True)
        y_space = np.linspace(self.down, self.up, num = 10, endpoint=True)
        inside = []
        for x in x_space:
            for y in y_space:
                p = Position(x,y)
#                print(p)
                if self.isInside(p):
                    inside.append(p)
        return inside
    
###############################################################################
    
    # -------------------------------------------------------------------------
    # -- CLASS Animation
    # -------------------------------------------------------------------------

class Animation:
    def __init__(self,
                 ballset = [],
                 gravity = 9.8,
                 space = [-2,2,2,-2]): #gauche, droite, haut, bas
        self.width = space[1] - space[0]
        self.height = space[2] - space[3]
        self.ballset = ballset
        self.gravity = gravity
        self.space = space
        self.time = 0
    
    def animate(self,dt,force = None):
        self.time += dt
        
        def zeroForce(x):
            return 0,0
        
        if force is None:
            force = zeroForce
        
        # update
        for ball in self.ballset:
            
            # speed
            forceX, forceY = force(ball.skeleton)
            forceY -= self.gravity * ball.getMass()
            ball.updateSpeed(ball.getVX() + forceX * dt / ball.getMass(),
                             ball.getVY() + forceY * dt / ball.getMass())
            
            # position
            new_x = ball.getX() + ball.getVX() * dt
            new_y = max(self.space[3]+0.75*ball.size, ball.getY() + ball.getVY() * dt)
            newPos = Position(new_x,new_y)
            ball.updatePos(newPos)
        
        # delect colliding pairs
        for i in range(len(self.ballset)):
            for j in range(i):
                if self.ballset[i].isColliding(self.ballset[j]):
                    
                    # update speed
                    vx1 = self.ballset[i].getVX()
                    vy1 = self.ballset[i].getVY()
                    m1 = self.ballset[i].getMass()
                    vx2 = self.ballset[j].getVX()                 
                    vy2 = self.ballset[j].getVY()
                    m2 = self.ballset[j].getMass()
                    
                    ux1 = (vx1 * (m1 - m2) + 2 * m2 * vx2)/(m1 + m2)
                    uy1 = (vy1 * (m1 - m2) + 2 * m2 * vy2)/(m1 + m2)
                    ux2 = (vx2 * (m2 - m1) + 2 * m1 * vx1)/(m1 + m2)
                    uy2 = (vy2 * (m2 - m1) + 2 * m1 * vy1)/(m1 + m2)
                    
                    self.ballset[i].updateSpeed(ux1,uy1)
                    self.ballset[j].updateSpeed(ux2,uy2)
         
        # detect collision with the walls
        for ball in self.ballset:
            x, y = ball.isBouncing(self.space)
            if x:
                ball.bounceX()
            if y:
                ball.bounceY()
    
    def _setdata(self):
        # determine the points to be diplayed
        self.ballpoints = []
        for k in range(len(self.ballset)):
            ball = self.ballset[k]
            bp = ball.getInside()
            self.ballpoints.append((bp,ball.color))
                        
                    
            

# init
n = 30
scale_speed = 10
border_size = 20
dt = 1. / 30 # 30fps

aux = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
colors = []
for col in aux.values():
    colors.append(col)
nb_colors = len(colors)

random.seed(0)
bs = []
for k in range(n):
    x = random.random() - 0.5
    y = random.random() - 0.5
    vx = (random.random() - 0.5) * scale_speed
    vy = (random.random() - 0.5) * scale_speed
    p = Position(x=x,y=y)
    s = Skeleton(position=p, vx=vx, vy=vy)
    b = Ball(skeleton = s, border_size = border_size,color = colors[k % nb_colors])
    bs.append(b)

animation = Animation(bs)

# set up the plotting values
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-3.2, 3.2), ylim=(-2.4, 2.4))

# determine the container in which the balls are bouncing
point = [animation.space[0], animation.space[3]]
container = plt.Rectangle(point, animation.width, animation.height,
                     ec='none', lw=2, fc='none')
ax.add_patch(container)
ms = fig.get_figwidth() / 200

balls_aux = []
for k in range(n):
    b, = ax.plot([],[],color = colors[k % nb_colors],ms=1)
    balls_aux.append(b)
balls = tuple(balls_aux)

def _init_function(): # used to initialize the animation
    global animation, container, balls, ball1, ball2, ball3, ball4, ball5
    
    for b in balls:
        b.set_data([],[])
#    ball1.set_data([],[])
#    ball2.set_data([],[])
#    ball3.set_data([],[])
#    ball4.set_data([],[])
#    ball5.set_data([],[])
    
    container.set_edgecolor('k')
#    return ball1, ball2, ball3, ball4, ball5, container
    return balls

def _animate(t): # one step of animation
    global animation, container, ball1, ball2, ball3, ball4, ball5, dt, fig, ax, ms
    animation.animate(dt)
    animation._setdata()
    
    for k in range(len(balls)):
        pos,col = animation.ballpoints[k]
        x_pos = [p.x for p in pos]
        y_pos = [p.y for p in pos]
        balls[k].set_data(x_pos,y_pos)
        balls[k].set_markersize(0.05)
#    
#    pos,col = animation.ballpoints[0]
##    print("ball1 : nb vert = {}".format(len(pos)))
#    x_pos = [p.x for p in pos]
#    y_pos = [p.y for p in pos]
#    ball1.set_data(x_pos,y_pos)
#    ball1.set_markersize(0.05)
#    
#    pos,col = animation.ballpoints[1]
##    print("ball2 : nb vert = {}".format(len(pos)))
#    x_pos = [p.x for p in pos]
#    y_pos = [p.y for p in pos]
#    ball2.set_data(x_pos,y_pos)
#    ball2.set_markersize(0.05)
#    
#    pos,col = animation.ballpoints[2]
##    print("ball3 : nb vert = {}".format(len(pos)))
#    x_pos = [p.x for p in pos]
#    y_pos = [p.y for p in pos]
#    ball3.set_data(x_pos,y_pos)
#    ball3.set_markersize(0.05)
#    
#    pos,col = animation.ballpoints[3]
##    print("ball4 : nb vert = {}".format(len(pos)))
#    x_pos = [p.x for p in pos]
#    y_pos = [p.y for p in pos]
#    ball4.set_data(x_pos,y_pos)
#    ball4.set_markersize(0.05)
#    
#    pos,col = animation.ballpoints[4]
##    print("ball5 : nb vert = {}".format(len(pos)))
#    x_pos = [p.x for p in pos]
#    y_pos = [p.y for p in pos]
#    ball5.set_data(x_pos,y_pos)
#    ball5.set_markersize(0.05)
#
#    return ball1, ball2, ball3, ball4, ball5, container
    return balls

res = anim.FuncAnimation(fig, _animate, frames=1000,
                              interval=10, blit=True, init_func=_init_function)

res.save('movement.mp4', fps=20, extra_args=['-vcodec', 'libx264'])

plt.show()

print("done")