import numpy
import cvxopt

class Cvxopt_Solve:
    def __init__(self, marg_xy, marg_xz):
        # marg_xy is a dictionary     (x,y) --> positive double
        # marg_xz is a dictionary     (x,z) --> positive double

        self.orig_marg_xy = None
        self.orig_marg_xz = None
        self.q_xy         = None
        self.q_xz         = None
        self.var_idx      = None
        self.X            = None
        self.Y            = None
        self.Z            = None
        self.A            = None
        self.b            = None
        self.marg_of_idx  = None  # triples xyz -- the missing one is None
        self.create_equations_called = False
        self.G            = None
        self.h            = None
        self.create_ieqs_called = False
        self.p_0          = None # initial solution
        self.solver_ret   = None
        self.p_final      = None
        self.set_to_zero  = _set_to_zero

        # Actual code:
        self.orig_marg_xy = dict(marg_xy)
        self.orig_marg_xz = dict(marg_xz)
        self.X = set( [ x   for x,y in self.orig_marg_xy.keys() ] + [ x   for x,z in self.orig_marg_xz.keys() ] )
        self.Y = set( [  y  for x,y in self.orig_marg_xy.keys() ] )
        self.Z = set(                                               [  z  for x,z in self.orig_marg_xz.keys() ] )
    # __init__()

    # tidy_up_distrib():
    def tidy_up_distrib(self,p):
        eps = 0
        # returns a tidied-up copy of p: every entry smaller than eps is treated as 0.
        p_new = dict()
        one = 0.
        for x,r in p.items():
            if r>=eps:
                p_new[x] = r
                one += r
        # Re-normalize --- Should I drop this ??? ?!?
        for x,r in p_new.items():
            p_new[x] = r / one;

        return p_new
    #^ tidy_up_distrib()

    # create_equations():
    def create_equations(self):
        # The point is that all the double values are considered non-zero
        # This function
        #  - creates the sets X,Y,Z
        #  - creates the dictionary var_idx: variables --> index of the variable
        #  - creates the matrices A,b: two types of marginal equations: p(x,y,*)==p_xy(x,y), p(x,*,z)==p_xz(x,z)
        #  - A has full rank: unneeded eqns are thrown out
        if self.create_equations_called:
            print("Some dork called create_equations twice...")
            exit(1)
        self.create_equations_called = True

        count_vars = 0
        self.var_idx    = dict()
        for x in self.X:
            for y in self.Y:
                if (x,y) in self.q_xy.keys():
                    for z in self.Z:
                        if (x,y,z) not in self.set_to_zero  and  (x,z) in self.q_xz.keys():
                            self.var_idx[ (x,y,z) ] = count_vars
                            count_vars += 1

        list_b           = [] # list of RHSs
        list_At          = [] # coefficient matrix in row-major (need column-major later)
        list_At_throwout = [] # DEBUG: omited equations go here---then we check whether the ranks are equal
        numo_thrownout   = 0

        self.marg_of_idx = []
        # xy-marginal equations
        for xy,rhs in self.q_xy.items():
            x,y = xy
            a = [ 0   for xyz in self.var_idx.keys() ] # initialize the whole row with 0
            for z in self.Z:
                if (x,y,z) in self.var_idx.keys():
                    i  = self.var_idx[ (x,y,z) ]
                    a[i] = 1.
            # list_At += a # splice !
            # list_b.append( rhs )
            # self.marg_of_idx.append(  (x,y,None)  )

            # We test if adding this equation increases the rank.
            # Because of the deleted variables (in set_to_zero), I don't know of a better way to do this...
            tmp_At = matrix( list_At+a, ( len(self.var_idx),  len(list_b)+1 ), 'd' )
            if numpy.linalg.matrix_rank( tmp_At ) > len(list_b):
                list_At += a # splice !
                list_b.append( rhs )
                self.marg_of_idx.append(  (x,None,z)  )

        # xz-marginal equations
        for xz,rhs in self.q_xz.items():
            x,z = xz
            a = [ 0   for xyz in self.var_idx.keys() ] # initialize the whole row with 0
            for y in self.Y:
                if (x,y,z) in self.var_idx.keys():
                    i = self.var_idx[ (x,y,z) ]
                    a[i] = 1.
            # Rank-check again:
            tmp_At = matrix( list_At+a, ( len(self.var_idx),  len(list_b)+1 ), 'd' )
            if numpy.linalg.matrix_rank( tmp_At ) > len(list_b):
                list_At += a # splice !
                list_b.append( rhs )
                self.marg_of_idx.append(  (x,None,z)  )

        # Now we create the CvxOpt matrix.
        self.b  = matrix( list_b,                   ( len(list_b),        1                          ), 'd' )
        At      = matrix( list_At,                  ( len(self.var_idx),  len(list_b)                ), 'd' )
        self.A = At.T
        rk = numpy.linalg.matrix_rank(self.A)
        if ( rk != len(list_b) ):
            print("BUG: There's something wrong with the rank of the coefficient matrix: it is ",rk," it should be ",len(list_b))
            exit(1)
        dim_space = len(self.var_idx)-rk
        print("Solution space has dimension ",dim_space)
    #^ create_equations()

    # create_ieqs():
    def create_ieqs(self):
        if not self.create_equations_called:
            print("You have to call create_equations() before calling create_ieqs()")
            exit(1)
        if self.create_ieqs_called:
            print("Some dork called create_ieqs() twice...")
            exit(1)
        self.create_ieqs_called = True
        self.G = spdiag( matrix( -1., (len(self.var_idx),1), 'd' ) )
        self.h = matrix( 0, (len(self.var_idx),1), 'd' )
    #^ create_ieqs()

    # CALLBACK fn for computing  f, grad f, Hess f
    def callback(self,p=None, zz=None):
        N = len(self.var_idx)
        if p is None:
            list_p_0 = [ 0.   for xyz in self.var_idx.keys() ]
            for xyz,i in self.var_idx.items():
                list_p_0[i] = 1. # self.p_0[xyz]
            # This is returns the starting solution for the iterative solution of the CP --- this is the 1st point to experiment with other distributions with the same marginals.
            return 0, matrix(list_p_0, (N,1), 'd' )

        # check if p is in the feasible region for the objective function:
        if min(p) <= 0 or max(p) > 1:
            return None

        p_dict = dict( (xyz,p[i]) for xyz,i in self.var_idx.items() )
        p_yz = marginal_yz(p_dict)

        # Compute f(p)
        f = 0
        for xyz,i in self.var_idx.items():
            x,y,z = xyz
            if p[i] > 0: f += p[i]*log(p[i]/p_yz[y,z])
#            if p[i] > 1.e-10: f += p[i]*log(p[i]/p_yz[y,z])

        # Compute gradient-transpose Df(p)
        list_Df = [ 0. for xyz in self.var_idx.keys() ]
        for xyz,i in self.var_idx.items():
            x,y,z = xyz
            if p[i] > 0:  list_Df[i] = log( p[i] / p_yz[y,z] )
#            if p[i] > 1.e-30:  list_Df[i] = log( p[i] / p_yz[y,z] )
        Df = matrix(list_Df, (1,N), 'd')

        if zz is None:
            return f,Df

        # Compute zz[0] * Hess f
        # This will be a sparse matrix
        entries = []
        rows    = []
        columns = []
        for xyz,i in self.var_idx.items():
            x,y,z = xyz
            p_yz__x = p_yz[y,z] - p[i] # sum_{* \ne x} p(*,y,z).
            for x_ in self.X:
                if x_==x: # diagonal
                    rows.append( i )
                    columns.append( i )
                    tmp_quot = zz[0] * p_yz__x / p_yz[y,z] # 1/p[x,y,z] - 1/p[*,y,z] = ( p[*,y,z] - p[x,y,z] )/( p[*,y,z] p[x,y,z] )
                    if p[i] > 0 * tmp_quot:   entries.append( tmp_quot / p[i] )
#                    if p[i] > 1.e-100 * tmp_quot:   entries.append( tmp_quot / p[i] )
                    else:
                        print("TROUBLE computing Hessian (diagonal)")
                        entries.append( 0. )
                else: # off diagonal
                    if (x_,y,z) in self.var_idx:
                        j = self.var_idx[ (x_,y,z) ]
                        val = - zz[0] / p_yz[y,z]
                        rows.append( i )
                        columns.append( j )
                        entries.append( val )
                # if diagonal
            # for x_
        # for xyz,i

        zH = spmatrix( entries, rows, columns, (N,N), 'd')
        if self.verbose_output: print("p=",list(p))
        return f,Df,zH
    #^ callback()


    # make_initial_solution()
    def make_initial_solution(self):
        for xyz in self.var_idx.keys():
            self.p_0[xyz] = 1.
    #^ make_initial_solution()

    # solve_it():
    def solve_it(self):
        self.q_xy = self.tidy_up_distrib(self.orig_marg_xy)
        self.q_xz = self.tidy_up_distrib(self.orig_marg_xz)

        self.create_equations()
        self.create_ieqs()
        self.make_initial_solution()

        self.solver_ret   = solvers.cp(self.callback, G=self.G, h=self.h, A=self.A, b=self.b)
        print("Solver terminated with status ",self.solver_ret['status'])
    #^ solve_it()

    def check_feasibility(self):
        status = self.solver_ret['status'] # This is a string

        # Make q
        q = dict()
        for xyz,i in self.var_idx.items():
            p[xyz] = self.solver_ret['x'][i]
        self.p_final = q

        obj_val,gradient = self.callback(p)

        # nonnegativity

        self.A*p - self.b

        # marginals_1
        # marginals_2
        # marginals_Inf

        llambda = self.solver_ret['z'] # array; fits with A I guess...

        mu = grad - self.A.T * llambda
        # mu_nonneg_viol

        # complementarity_max
        # complementarity_sum

        return ..., ..., ..., ..., ....
    #^ check_feasibility()

#^ class Cvxopt_Solve
