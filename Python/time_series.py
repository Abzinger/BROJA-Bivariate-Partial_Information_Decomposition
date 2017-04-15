#!/usr/bin/python3
import TartuSynergy
import sys
import random
import json
#
# Main function: main()
# =====================
#
# Usage:
#             main([0, fun, noise_level, n_samples])
#        or   main([0, '-v', fun, noise_level, n_samples])
#
# Makes n_samples samples from function `fun' (fun is a string
# identifying the function), adds noise to it, then determines shared
# information and synergetic information.
# If -v is given, write lots and lots of progress output.
#
# Add new functions
# =================
#
# To add a new function, you have to define a class, and add an object
# of it to the dictionary of available functions.  (This dictionary
# maps capitalized names of function to function objects.)
#
# The functions are defined throuh classes.
# The class must containt the following
#
# *  n_u, n_v, n_w --- desired input length
#
# *  fun(u,v,w) -- take three independent random variables
#    u in {0,...,n_u-1}, v in {0,...,n_v-1}, w in {0,...,n_w-1}
#    Return a triple X,Y,Z
#
# *  noise(X) adds random noise to the X-part of X,Y,Z
#
# * true_input_distrib()    true PDF (not sampled)
#   true_result_distrib()   right optimal PDF (when input is not sampled)
#   true_CI()               true synergetic info (when input is not sampled)
#   true_SI()               true shared info  (when input is not sampled)


#^ END OF DOCUMENTATION

the_available_functions = dict()

class Fun_Data_rdn:
    def __init__(self):
        self.n_u=self.n_v=0
        self.n_w=2
    def fun(self,y,z,w):
        return w,w,w
    def noise(self,x):
        x  = ( x +random.randint(0,1) )%2
        return x
    def true_input_distrib(self):
        return {(0,0,0):.5, (1,1,1):.5}
    def true_result_distrib(self):
        return {(0,0,0):.5, (1,1,1):.5}
    def true_CI(self):
        return 0
    def true_SI(self):
        return 1


the_available_functions["RDN"] = Fun_Data_rdn()


class Fun_Data_unq:
    def __init__(self):
        self.n_u=self.n_v=2
        self.n_w=0
    def fun(self,y,z,w):
        return (y,z),y,z
    def noise(self,x):
        y,z  = x
        y = ( y +random.randint(0,1) )%2
        z = ( z +random.randint(0,1) )%2
        return y,z
    def true_input_distrib(self):
        return {((0,0),0,0):.25, ((0,1),0,1):.25, ((1,0),1,0):.25, ((1,1),1,1):.25}
    def true_result_distrib(self):
        return {((0,0),0,0):.25, ((0,1),0,1):.25, ((1,0),1,0):.25, ((1,1),1,1):.25}
    def true_CI(self):
        return 0 # this is a mistake in BROJA'14
    def true_SI(self):
        return 0
the_available_functions["UNQ"] = Fun_Data_unq()


class Fun_Data_xor:
    def __init__(self):
        self.n_u=self.n_v=2
        self.n_w=0
    def fun(self,y,z,w):
        return (y+z)%2,y,z
    def noise(self,x):
        x  = ( x +random.randint(0,1) )%2
        return x
    def true_input_distrib(self):
        return {(0,0,0):.25, (1,0,1):.25, (1,1,0):.25, (0,1,1):.25}
    def true_result_distrib(self):
        return {(0,0,0):.25, (1,0,1):.25, (1,1,0):.25, (0,1,1):.25}
    def true_CI(self):
        return 1
    def true_SI(self):
        return 0
the_available_functions["XOR"] = Fun_Data_xor()



class Fun_Data_and:
    def __init__(self):
        self.n_u=self.n_v=2
        self.n_w=0
    def fun(self,y,z,w):
        return (y*z)%2,y,z
    def noise(self,x):
        x  = ( x +random.randint(0,1) )%2
        return x
    def true_input_distrib(self):
        return {(0,0,0):.25, (0,0,1):.25, (0,1,0):.25, (1,1,1):.25}
    def true_result_distrib(self):
        return { (0,0,0): .5, (0,1,1): .25, (1,1,1): .25 }
    def true_CI(self):
        return .5
    def true_SI(self):
        return .311278124459132843017578125 # I hope
the_available_functions["AND"] = Fun_Data_and()


class Fun_Data_rdnxor:
    def __init__(self):
        self.n_u=self.n_v=self.n_w=2
    def fun(self,y,z,w):
        return ((y+z)%2,w), (y,w), (z,w)
    def noise(self,x):
        x1,w = x
        x1 = ( x1+random.randint(0,1) )%2
        w  = ( w +random.randint(0,1) )%2
        return x1,w
    def true_input_distrib(self):
        return {
            ((0,0), (0,0), (0,0)): .125,
            ((1,0), (0,0), (1,0)): .125,
            ((1,0), (1,0), (0,0)): .125,
            ((0,0), (1,0), (1,0)): .125,
            ((0,1), (0,1), (0,1)): .125,
            ((1,1), (0,1), (1,1)): .125,
            ((1,1), (1,1), (0,1)): .125,
            ((0,1), (1,1), (1,1)): .125
        }
    def true_result_distrib(self):
        return {((0, 0), (0, 0), (0, 0)):0.0625, ((0, 0), (0, 0), (1, 0)):0.0625, ((0, 0), (1, 0), (0, 0)):0.0625, ((0, 0), (1, 0), (1, 0)):0.0625, ((0, 1), (0, 1), (0, 1)):0.0625, ((0, 1), (0, 1), (1, 1)):0.0625, ((0, 1), (1, 1), (0, 1)):0.0625, ((0, 1), (1, 1), (1, 1)):0.0625, ((1, 0), (0, 0), (0, 0)):0.0625, ((1, 0), (0, 0), (1, 0)):0.0625, ((1, 0), (1, 0), (0, 0)):0.0625, ((1, 0), (1, 0), (1, 0)):0.0625, ((1, 1), (0, 1), (0, 1)):0.0625, ((1, 1), (0, 1), (1, 1)):0.0625, ((1, 1), (1, 1), (0, 1)):0.0625, ((1, 1), (1, 1), (1, 1)):0.0625}
    def true_CI(self):
        return 1
    def true_SI(self):
        return 1
the_available_functions["RDNXOR"] = Fun_Data_rdnxor()

class Fun_Data_rdnunqxor:
    def __init__(self):
        self.n_u=self.n_v=4
        self.n_w=2
    def fun(self,y,z,w):
        y1 = y%2
        y2 = y//2
        z1 = z%2
        z2 = z//2
        return ((y1+z1)%2,y2,z2,w), (y,w), (z,w)
    def noise(self,x):
        x1,y2,z2,w = x
        x1 = ( x1+random.randint(0,1) )%2
        y2 = ( y2+random.randint(0,1) )%2
        z2 = ( z2+random.randint(0,1) )%2
        w  = ( w +random.randint(0,1) )%2
        return x1,y2,z2,w
    def true_input_distrib(self):
        return {
            ((0,0,0,0), (0,0), (0,0)): 1/32,
            ((1,0,0,0), (0,0), (1,0)): 1/32,
            ((0,0,1,0), (0,0), (2,0)): 1/32,
            ((1,0,1,0), (0,0), (3,0)): 1/32,
            ((1,0,0,0), (1,0), (0,0)): 1/32,
            ((0,0,0,0), (1,0), (1,0)): 1/32,
            ((1,0,1,0), (1,0), (2,0)): 1/32,
            ((0,0,1,0), (1,0), (3,0)): 1/32,
            ((0,1,0,0), (2,0), (0,0)): 1/32,
            ((1,1,0,0), (2,0), (1,0)): 1/32,
            ((0,1,1,0), (2,0), (2,0)): 1/32,
            ((1,1,1,0), (2,0), (3,0)): 1/32,
            ((1,1,0,0), (3,0), (0,0)): 1/32,
            ((0,1,0,0), (3,0), (1,0)): 1/32,
            ((1,1,1,0), (3,0), (2,0)): 1/32,
            ((0,1,1,0), (3,0), (3,0)): 1/32,
            ((0,0,0,1), (0,1), (0,1)): 1/32,
            ((1,0,0,1), (0,1), (1,1)): 1/32,
            ((0,0,1,1), (0,1), (2,1)): 1/32,
            ((1,0,1,1), (0,1), (3,1)): 1/32,
            ((1,0,0,1), (1,1), (0,1)): 1/32,
            ((0,0,0,1), (1,1), (1,1)): 1/32,
            ((1,0,1,1), (1,1), (2,1)): 1/32,
            ((0,0,1,1), (1,1), (3,1)): 1/32,
            ((0,1,0,1), (2,1), (0,1)): 1/32,
            ((1,1,0,1), (2,1), (1,1)): 1/32,
            ((0,1,1,1), (2,1), (2,1)): 1/32,
            ((1,1,1,1), (2,1), (3,1)): 1/32,
            ((1,1,0,1), (3,1), (0,1)): 1/32,
            ((0,1,0,1), (3,1), (1,1)): 1/32,
            ((1,1,1,1), (3,1), (2,1)): 1/32,
            ((0,1,1,1), (3,1), (3,1)): 1/32
        }
    def true_result_distrib(self):
        return {((0, 0, 0, 0), (0, 0), (0, 0)):1/64, ((0, 0, 0, 0), (0, 0), (1, 0)):1/64, ((0, 0, 0, 0), (1, 0), (0, 0)):1/64, ((0, 0, 0, 0), (1, 0), (1, 0)):1/64, ((0, 0, 0, 1), (0, 1), (0, 1)):1/64, ((0, 0, 0, 1), (0, 1), (1, 1)):1/64, ((0, 0, 0, 1), (1, 1), (0, 1)):1/64, ((0, 0, 0, 1), (1, 1), (1, 1)):1/64, ((0, 0, 1, 0), (0, 0), (2, 0)):1/64, ((0, 0, 1, 0), (0, 0), (3, 0)):1/64, ((0, 0, 1, 0), (1, 0), (2, 0)):1/64, ((0, 0, 1, 0), (1, 0), (3, 0)):1/64, ((0, 0, 1, 1), (0, 1), (2, 1)):1/64, ((0, 0, 1, 1), (0, 1), (3, 1)):1/64, ((0, 0, 1, 1), (1, 1), (2, 1)):1/64, ((0, 0, 1, 1), (1, 1), (3, 1)):1/64, ((0, 1, 0, 0), (2, 0), (0, 0)):1/64, ((0, 1, 0, 0), (2, 0), (1, 0)):1/64, ((0, 1, 0, 0), (3, 0), (0, 0)):1/64, ((0, 1, 0, 0), (3, 0), (1, 0)):1/64, ((0, 1, 0, 1), (2, 1), (0, 1)):1/64, ((0, 1, 0, 1), (2, 1), (1, 1)):1/64, ((0, 1, 0, 1), (3, 1), (0, 1)):1/64, ((0, 1, 0, 1), (3, 1), (1, 1)):1/64, ((0, 1, 1, 0), (2, 0), (2, 0)):1/64, ((0, 1, 1, 0), (2, 0), (3, 0)):1/64, ((0, 1, 1, 0), (3, 0), (2, 0)):1/64, ((0, 1, 1, 0), (3, 0), (3, 0)):1/64, ((0, 1, 1, 1), (2, 1), (2, 1)):1/64, ((0, 1, 1, 1), (2, 1), (3, 1)):1/64, ((0, 1, 1, 1), (3, 1), (2, 1)):1/64, ((0, 1, 1, 1), (3, 1), (3, 1)):1/64, ((1, 0, 0, 0), (0, 0), (0, 0)):1/64, ((1, 0, 0, 0), (0, 0), (1, 0)):1/64, ((1, 0, 0, 0), (1, 0), (0, 0)):1/64, ((1, 0, 0, 0), (1, 0), (1, 0)):1/64, ((1, 0, 0, 1), (0, 1), (0, 1)):1/64, ((1, 0, 0, 1), (0, 1), (1, 1)):1/64, ((1, 0, 0, 1), (1, 1), (0, 1)):1/64, ((1, 0, 0, 1), (1, 1), (1, 1)):1/64, ((1, 0, 1, 0), (0, 0), (2, 0)):1/64, ((1, 0, 1, 0), (0, 0), (3, 0)):1/64, ((1, 0, 1, 0), (1, 0), (2, 0)):1/64, ((1, 0, 1, 0), (1, 0), (3, 0)):1/64, ((1, 0, 1, 1), (0, 1), (2, 1)):1/64, ((1, 0, 1, 1), (0, 1), (3, 1)):1/64, ((1, 0, 1, 1), (1, 1), (2, 1)):1/64, ((1, 0, 1, 1), (1, 1), (3, 1)):1/64, ((1, 1, 0, 0), (2, 0), (0, 0)):1/64, ((1, 1, 0, 0), (2, 0), (1, 0)):1/64, ((1, 1, 0, 0), (3, 0), (0, 0)):1/64, ((1, 1, 0, 0), (3, 0), (1, 0)):1/64, ((1, 1, 0, 1), (2, 1), (0, 1)):1/64, ((1, 1, 0, 1), (2, 1), (1, 1)):1/64, ((1, 1, 0, 1), (3, 1), (0, 1)):1/64, ((1, 1, 0, 1), (3, 1), (1, 1)):1/64, ((1, 1, 1, 0), (2, 0), (2, 0)):1/64, ((1, 1, 1, 0), (2, 0), (3, 0)):1/64, ((1, 1, 1, 0), (3, 0), (2, 0)):1/64, ((1, 1, 1, 0), (3, 0), (3, 0)):1/64, ((1, 1, 1, 1), (2, 1), (2, 1)):1/64, ((1, 1, 1, 1), (2, 1), (3, 1)):1/64, ((1, 1, 1, 1), (3, 1), (2, 1)):1/64, ((1, 1, 1, 1), (3, 1), (3, 1)):1/64}
    def true_CI(self):
        return 1
    def true_SI(self):
        return 1
the_available_functions["RDNUNQXOR"] = Fun_Data_rdnunqxor()


class Fun_Data_xorand:
    def __init__(self):
        self.n_u=self.n_v=2
        self.n_w=0
    def fun(self,y,z,w):
        return ((y+z)%2,y*z), y, z
    def noise(self,x):
        x1,x2 = x
        return ( x1+random.randint(0,1) )%2, ( x2+random.randint(0,1) )%2
    def true_input_distrib(self):
        return {
            ((0,0), 0, 0): .25,
            ((1,0), 0, 1): .25,
            ((1,0), 1, 0): .25,
            ((0,1), 1, 1): .25,
        }
    def true_result_distrib(self):
        return {((0, 0), 0, 0):0.25, ((0, 1), 1, 1):0.25, ((1, 0), 0, 0):0.25, ((1, 0), 1, 1):0.25}
    def true_CI(self):
        return 1
    def true_SI(self):
        return 1/2
the_available_functions["XORAND"] = Fun_Data_xorand()

def to_list(x):
    if x == ():
        return []
    if type(x) != tuple:
        return x
    return [to_list(x[0])] + to_list(x[1:])

def test__solve_time_series(filename, fun_obj, noise_level,   numo_samples, verbose=False):
    y_list=[]
    z_list=[]
    x_list=[]
    counts = dict()

    if numo_samples == 0:
        pdf = fun_obj.true_input_distrib()
    else:
        L = numo_samples

        for i in range(0, L):
            if fun_obj.n_u>0: new_y = random.randint(0,fun_obj.n_u-1)
            else:             new_y = None
            if fun_obj.n_v>0: new_z = random.randint(0,fun_obj.n_v-1)
            else:             new_z = None
            if fun_obj.n_w>0: new_w = random.randint(0,fun_obj.n_w-1)
            else:             new_w = None
            x,y,z = fun_obj.fun(new_y,new_z,new_w)
            x_list.append( x )
            y_list.append( y )
            z_list.append( z )

        # add noise
        for i in range(0,L):
            noise = 0
            p = random.random();
            if p<noise_level:
                x_list[i] = fun_obj.noise( x_list[i] )

        # update counts
        for i in range(0,L):
            if (x_list[i],y_list[i],z_list[i]) in counts.keys():   counts[(x_list[i],y_list[i],z_list[i])] += 1
            else:                                                  counts[(x_list[i],y_list[i],z_list[i])] =  1

        # make pdf from counts
        pdf = dict()
        for xyz,c in counts.items():  pdf[xyz] = c/float(L)
        # Finished sampling creation of pdf
    #END of creation of pdf

    if filename != None:
        with open(filename, 'w') as f:
            array = []
            for k,v in pdf.items():
                array.append(   [ to_list(k), v ]   )
            json.dump(array, f)
        return 0
    else:
        print("pdf=",TartuSynergy.sorted_pdf(pdf))
        print("test__solve_time_series(): L = ",numo_samples)
        retval=TartuSynergy.solve_PDF( pdf, fun_obj.true_input_distrib(), fun_obj.true_result_distrib(), fun_obj.true_CI(), fun_obj.true_SI() , verbose=verbose)
        optpdf,feas,kkt,CI,SI = retval
        print("CI:",CI,"  SI:",SI,"  sum of marginal eqns violations:",feas,"  maximal KKT-system constraint violation:",kkt)
        return retval
    #
#^ test__solve_time_series()



####################################################################################################

def print_parameters_usage_msg():
    print("Usage:   time_series.py -h")
    print("   or:   time_serie.py  [-v|w] fun noiselevel # ")
    print("time_series.py -h    lists the functions.")
    print("Otherwise,")
    print("        fun          is a function;")
    print("        noiselevel   is the probability of noise, in [0,1[;")
    print("        #            is number of samples; with 0 the true distribution is used.")
    print("Unless -w is given, the information decomposition is computed.")
    print("In that case:")
    print("        -v           sets verbose mode")
    print("With -w, the probability density function is written to")
    print("the file `fun-noiselevel-#.dens'")

def main(sys_argv):
    verbose_mode = False
    file_mode    = False
    if len(sys_argv)==5:
        if sys_argv[1]=='-v':
            verbose_mode = True
        elif sys_argv[1]=='-w':
            file_mode = True
        else:
            print("There's something wrong with the parameters!")
            print_parameters_usage_msg()
            return 1
        del sys_argv[1]

    if len(sys_argv) <= 1:
        print_parameters_usage_msg()
        return 0
    elif len(sys_argv) == 2:
        if sys_argv[1]=='-h':
            print("Available functions:")
            print(the_available_functions.keys())
            return 0
        else:
            print("There's something wrong with the parameters!")
            print_parameters_usage_msg()
            return 1
    elif len(sys_argv) != 4:
        print("There's something wrong with the parameters!")
        print_parameters_usage_msg()
        return 1
    else:
        fun_S        = (str(sys_argv[1])).upper()
        n_u=n_v=n_w  = None
        noise_p      = float(sys_argv[2])
        numo_samples = int(sys_argv[3])
        fun_obj      = None
        true_input  = None
        true_result = None
        true_CI     = None
        true_SI     = None
        valid_parameters  = True
        if fun_S in the_available_functions.keys():
            fun_obj = the_available_functions[fun_S]
        else:
            print("I don't recognize the function ", fun_S)
            print("Available functions:")
            print(the_available_functions.keys())
            valid_parameters = False

        if noise_p<0 or noise_p >= 1:
            print("Something is wrong with the noise probability ", noise_p)
            print_parameters_usage_msg()
            valid_parameters = False
        if numo_samples<0 or numo_samples > 1.e20:
            print("Something is wrong with the number of samples", numo_samples)
            print_parameters_usage_msg()
            valid_parameters = False
        if valid_parameters:
#           solvers.options['show_progress'] = False
            if file_mode:
                filename = fun_S+"-"+str(noise_p)+"-"+str(numo_samples)+".dens"
            else:
                filename = None
            retval = test__solve_time_series(filename, fun_obj, noise_p, numo_samples, verbose=verbose_mode)
            return retval
        else:
            return 1
    #else main work
#^ main()

main(sys.argv)
