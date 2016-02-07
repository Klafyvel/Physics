import numpy as np
import matplotlib.pyplot as pl

GRAVITY = 9.81
LAMBDA = 0.5

class Part:
    def __init__(self, function=lambda x:0, param='CARTESIAN', xlim=(0,0), origin=(0,0)):
        """
        Initialize the part of an object.
        :param function: The function tha desribes the part
        :param param: 'CARTESIAN' | 'POLAR'
        :param xlim: The domain of definition.
        :param polar_origin: The origin of param.
        """

        self.function = function
        self.param = param
        self.xlim = xlim
        self.values = [[],[]]
        self.origin = origin
        self.last_test_point = None
        self.last_n_vect = None
        self.last_t_vect = None

    def get_points(self, res):
        """
        Computes the values of the function
        :param res: nb of points
        :return: list of points
        """
        if len(self.values[0]) is res:
            return self.values
        if self.param == 'CARTESIAN':
            self.values[0] = [x+self.origin[0] for x in np.linspace(self.xlim[0], self.xlim[1], res)]
            self.values[1] = [self.function(x)+self.origin[1] for x in np.linspace(self.xlim[0], self.xlim[1], res)]
        else:
            self.values = [[],[]]
            tetas = np.linspace(self.xlim[0], self.xlim[1], res)
            for t in tetas:
                r = self.function(t)
                self.values[0].append(np.cos(t)*r+self.origin[0])
                self.values[1].append(np.sin(t)*r+self.origin[1])
        return self.values
    def is_on_part_cartesian(self,point,precision):
        """
        Check if a point is on the part in cartesian mode.
        :param point: Your point
        :param precision: the scalar precision around the y value
        :return: None if not on part, the x value else
        """
        x,y = point
        if x < self.xlim[0] or x > self.xlim[1]:
            return False
        y_f = self.function(x-self.origin[0])+self.origin[1]
        #self.last_dist = abs(y-y_f)
        return x, abs(y-y_f)

    def is_on_part_polar(self, point, precision):
        """
        Check if a point is on the part in polar mode.
        :param point: Your point
        :param precision: the scalar precision around the y value
        :return: None if not on part, the tetha value else
        """
        x,y = point
        t_min,t_max = self.xlim
        m = (t_min+t_max)/2
        r=0
        while t_max-t_min > precision:
            r_min = self.function(t_min)
            r_max = self.function(t_max)
            x_f_max,y_f_max = r_max*np.cos(t_max)+self.origin[0],r_max*np.sin(t_max)+self.origin[1]
            x_f_min,y_f_min = r_min*np.cos(t_min)+self.origin[0],r_min*np.sin(t_min)+self.origin[1]
            d_max=((x_f_max-x)**2+(y_f_max-y)**2)**(1/2)
            d_min=((x_f_min-x)**2+(y_f_min-y)**2)**(1/2)

            if d_max < d_min:
                t_min = m
            else:
                t_max = m
            m = (t_max+t_min)/2
        r = self.function(m)
        x_f,y_f = r*np.cos(m)+self.origin[0],r*np.sin(m)+self.origin[1]
        d = ((x_f-x)**2+(y_f-y)**2)**1/2
        #self.last_dist = d
        return m, d

    def is_on_part(self,point,precision):
        """
        Check if a point is on the part.
        :param point: Your point
        :param precision: the scalar precision around the y value
        :return: None if not on part, the param value else
        """
        if self.param == 'CARTESIAN':
            return self.is_on_part_cartesian(point,precision)
        else:
            return self.is_on_part_polar(point,precision)


    def get_tangent_vector(self,point,displacement,precision=1):
        if precision ==0:
            precision = 10**(-15)
        #print(precision)
        inter = self.is_on_part(point, precision)
        if not inter:
            return (0,0)
        x,d = inter
        if d > precision:
            return (0,0)
        if self.param == 'CARTESIAN':
            self.last_test_point = (x,self.function(x)+self.origin[1])
            dx = 2*precision
            dy = self.function(x-self.origin[0]+precision)-self.function(x-self.origin[0]-precision)
        else:
            r1 = self.function(x-precision)
            r2 = self.function(x+precision)
            r = self.function(x)
            dx = np.cos(x+precision)*r2-np.cos(x-precision)*r1
            dy = np.sin(x+precision)*r2-np.sin(x-precision)*r1
            self.last_test_point = (np.cos(x)*r+self.origin[0], np.sin(x)*r+self.origin[1])
        d = (dx**2 + dy**2)**(1/2)
        dx,dy = dx/d, dy/d
        return (dx,dy)
    def get_resulting_force(self, pos, acting_force,velocity, precision, draw_target=None):
        t = self.get_tangent_vector(pos, velocity,precision)
        n = (-t[1], t[0])
        #if t[1]>0:
        #    print('bli')
        #    t = (t[0],-t[1])
        self.last_t_vect = t
        self.last_n_vect = n

        if draw_target and self.last_test_point:
            draw_target.arrow(self.last_test_point[0], self.last_test_point[1], t[0], t[1],head_width=1, color='orange')
            draw_target.arrow(self.last_test_point[0], self.last_test_point[1], n[0], n[1],head_width=1, color='orange')


        p_n = round(acting_force[0]*n[0] + acting_force[1]*n[1],5)

        r = (round(-n[0]*p_n,5),round(-n[1]*p_n,5))
        return r

class Point:
    def __init__(self,pos,mass,velocity=(0,0)):
        self.pos = pos
        self.mass = mass
        self.v = velocity
        self.a = (0,0)
        #self.draw = {}

        self.new_pos = None
        self.new_v = None

        self.forces = {}

    def move(self,objects,dt,draw_debug=None):
        """
        Moves the point.
        :param objects: The objects in the world
        :param dt: delta-time
        """
        p = (0,round(-self.mass * GRAVITY,5))
        r = (0,0)
        for o in objects:
            r_x,r_y = o.get_resulting_force(self.pos, p, self.v,precision=(self.v[0]**2+self.v[1]**2)**(1/2)*dt, draw_target=draw_debug)
            r = (round(r[0]+r_x,5),round(r[1]+r_y,5))

        print((p[0]+r[0],p[1]+r[1]))

        if (p[0]+r[0],p[1]+r[1])!=(0,0):
            dv = (round(dt / self.mass * (p[0]+r[0]),5),round(dt / self.mass * (p[1]+r[1]),5))
            self.new_v = (round(self.v[0]+dv[0],5), round(self.v[1]+dv[1],5))
        else:
            self.new_v = self.v[:]
        self.a = ((p[0]+r[0])/self.mass,(p[1]+r[1])/self.mass)
        dOM = (round(self.v[0]*dt,5), round(self.v[1]*dt,5))
        self.new_pos = (round(self.pos[0]+dOM[0],5), round(self.pos[1]+dOM[1],5))


        self.rectify(objects)


        self.forces = {
            '$\\vec{P}$':p,
            '$\\vec{R}$':r
        }

        """
        self.draw = {
            'M':self.pos,
            'force':{
                '$\\vec{P}$':p,
            },
            'velocity':self.v,
            'acceleration': (dv[0]/dt, dv[1]/dt)
        }
        if r != (0,0):
            self.draw['force']['$\\vec{R}$'] = r
        """

    def rectify(self,objects):
        #print('\n')
        if not self.new_pos or not self.new_v:
            return
        d_min = None
        i_d_min = None
        for i,o in enumerate(objects):
            if o.last_test_point:
                current_d = (self.pos[0]-o.last_test_point[0], self.pos[1]-o.last_test_point[1])
                d = (self.new_pos[0]-o.last_test_point[0], self.new_pos[1]-o.last_test_point[1])
                if round((o.last_n_vect[0]*d[0] + o.last_n_vect[1]*d[1])*(o.last_n_vect[0]*current_d[0] + o.last_n_vect[1]*current_d[1]),5) < 0:
                    #print("Dépassement sur",i)
                    if d_min:
                        d_min = min(d_min,(current_d[0]**2+current_d[1]**2)**(1/2))
                    else:
                        d_min = (current_d[0]**2+current_d[1]**2)**(1/2)
                    i_d_min = i
                    #print(i_d_min,d_min)


        if i_d_min:
            o = objects[i_d_min]
            print("Dépasse")
            self.pos = o.last_test_point[:]
            v = (self.new_v[0]+self.v[0])/2*o.last_t_vect[0] + (self.new_v[1]+self.v[1])/2*o.last_t_vect[1]
            self.v = (round(v*o.last_t_vect[0],5), round(v*o.last_t_vect[1],5))

        else:
            self.pos = self.new_pos[:]
            self.v = (round(self.new_v[0],5), round(self.new_v[1],5))
        self.new_pos, self.new_v = None, None

    def draw(self,subplot,color='#000000'):
        subplot.scatter(self.pos[0],self.pos[1],10,c='blue')
        for f in self.forces:
            v = self.forces[f]
            subplot.arrow(self.pos[0], self.pos[1], v[0]/3, v[1]/3, head_width=1, color=color)
        subplot.arrow(self.pos[0], self.pos[1], self.v[0], self.v[1],head_width=1, color='red')
        #subplot.arrow(self.pos[0], self.pos[1], self.a[0], self.a[1],head_width=1, color='green')


