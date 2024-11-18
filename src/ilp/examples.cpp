/* Copyright 2024, Gurobi Optimization, LLC */

/* Solve a traveling salesman problem on a randomly generated set of
   points using lazy constraints.   The base MIP model only includes
   'degree-2' constraints, requiring each node to have exactly
   two incident edges.  Solutions to this model may contain subtours -
   tours that don't visit every node.  The lazy constraint callback
   adds new constraints to cut them off. */

#define GUROBI_DIR "/opt/gurobi1103/linux64/include/gurobi_c++.h"

#include GUROBI_DIR
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>
using namespace std;

string itos(int i) {stringstream s; s << i; return s.str(); }
double distance(double* x, double* y, int i, int j);
void findsubtour(int n, double** sol, int* tourlenP, int* tour);

// Subtour elimination callback.  Whenever a feasible solution is found,
// find the smallest subtour, and add a subtour elimination constraint
// if the tour doesn't visit every node.

class subtourelim: public GRBCallback
{
  public:
    GRBVar** vars;
    int n;
    subtourelim(GRBVar** xvars, int xn) {
      vars = xvars;
      n    = xn;
    }
  protected:
    void callback() {
      try {
        if (where == GRB_CB_MIPSOL) {
          // Found an integer feasible solution - does it visit every node?
          double **x = new double*[n];
          int *tour = new int[n];
          int i, j, len;
          for (i = 0; i < n; i++)
            x[i] = getSolution(vars[i], n);

          findsubtour(n, x, &len, tour);

          if (len < n) {
            // Add subtour elimination constraint
            GRBLinExpr expr = 0;
            for (i = 0; i < len; i++)
              for (j = i+1; j < len; j++)
                expr += vars[tour[i]][tour[j]];
            addLazy(expr <= len-1);
          }

          for (i = 0; i < n; i++)
            delete[] x[i];
          delete[] x;
          delete[] tour;
        }
      } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }
    }
};

// Given an integer-feasible solution 'sol', find the smallest
// sub-tour.  Result is returned in 'tour', and length is
// returned in 'tourlenP'.

void
findsubtour(int      n,
            double** sol,
            int*     tourlenP,
            int*     tour)
{
  bool* seen = new bool[n];
  int bestind, bestlen;
  int i, node, len, start;

  for (i = 0; i < n; i++)
    seen[i] = false;

  start = 0;
  bestlen = n+1;
  bestind = -1;
  node = 0;
  while (start < n) {
    for (node = 0; node < n; node++)
      if (!seen[node])
        break;
    if (node == n)
      break;
    for (len = 0; len < n; len++) {
      tour[start+len] = node;
      seen[node] = true;
      for (i = 0; i < n; i++) {
        if (sol[node][i] > 0.5 && !seen[i]) {
          node = i;
          break;
        }
      }
      if (i == n) {
        len++;
        if (len < bestlen) {
          bestlen = len;
          bestind = start;
        }
        start += len;
        break;
      }
    }
  }

  for (i = 0; i < bestlen; i++)
    tour[i] = tour[bestind+i];
  *tourlenP = bestlen;

  delete[] seen;
}

// Euclidean distance between points 'i' and 'j'.

double
distance(double* x,
         double* y,
         int     i,
         int     j)
{
  double dx = x[i]-x[j];
  double dy = y[i]-y[j];

  return sqrt(dx*dx+dy*dy);
}

int
main(int   argc,
     char *argv[])
{
  if (argc < 2) {
    cout << "Usage: tsp_c++ size" << endl;
    return 1;
  }

  int n = atoi(argv[1]);
  double* x = new double[n];
  double* y = new double[n];

  int i;
  for (i = 0; i < n; i++) {
    x[i] = ((double) rand())/RAND_MAX;
    y[i] = ((double) rand())/RAND_MAX;
  }

  GRBEnv *env = NULL;
  GRBVar **vars = NULL;

  vars = new GRBVar*[n];
  for (i = 0; i < n; i++)
    vars[i] = new GRBVar[n];

  try {
    int j;

    env = new GRBEnv();
    GRBModel model = GRBModel(*env);

    // Must set LazyConstraints parameter when using lazy constraints

    model.set(GRB_IntParam_LazyConstraints, 1);

    // Create binary decision variables

    for (i = 0; i < n; i++) {
      for (j = 0; j <= i; j++) {
        vars[i][j] = model.addVar(0.0, 1.0, distance(x, y, i, j),
                                  GRB_BINARY, "x_"+itos(i)+"_"+itos(j));
        vars[j][i] = vars[i][j];
      }
    }

    // Degree-2 constraints

    for (i = 0; i < n; i++) {
      GRBLinExpr expr = 0;
      for (j = 0; j < n; j++)
        expr += vars[i][j];
      model.addConstr(expr == 2, "deg2_"+itos(i));
    }

    // Forbid edge from node back to itself

    for (i = 0; i < n; i++)
      vars[i][i].set(GRB_DoubleAttr_UB, 0);

    // Set callback function

    subtourelim cb = subtourelim(vars, n);
    model.setCallback(&cb);

    // Optimize model

    model.optimize();

    // Extract solution

    if (model.get(GRB_IntAttr_SolCount) > 0) {
      double **sol = new double*[n];
      for (i = 0; i < n; i++)
        sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);

      int* tour = new int[n];
      int len;

      findsubtour(n, sol, &len, tour);
      assert(len == n);

      cout << "Tour: ";
      for (i = 0; i < len; i++)
        cout << tour[i] << " ";
      cout << endl;

      for (i = 0; i < n; i++)
        delete[] sol[i];
      delete[] sol;
      delete[] tour;
    }

  } catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Error during optimization" << endl;
  }

  for (i = 0; i < n; i++)
    delete[] vars[i];
  delete[] vars;
  delete[] x;
  delete[] y;
  delete env;
  return 0;
}

/* Copyright 2024, Gurobi Optimization, LLC */
/*
  Sudoku example.

  The Sudoku board is a 9x9 grid, which is further divided into a 3x3 grid
  of 3x3 grids.  Each cell in the grid must take a value from 0 to 9.
  No two grid cells in the same row, column, or 3x3 subgrid may take the
  same value.

  In the MIP formulation, binary variables x[i,j,v] indicate whether
  cell <i,j> takes value 'v'.  The constraints are as follows:
    1. Each cell must take exactly one value (sum_v x[i,j,v] = 1)
    2. Each value is used exactly once per row (sum_i x[i,j,v] = 1)
    3. Each value is used exactly once per column (sum_j x[i,j,v] = 1)
    4. Each value is used exactly once per 3x3 subgrid (sum_grid x[i,j,v] = 1)

  Input datasets for this example can be found in examples/data/sudoku*.
*/

#include <sstream>
using namespace std;

#define sd 3
#define n (sd*sd)

string itos(int i) {stringstream s; s << i; return s.str(); }

int
main(int   argc,
     char *argv[])
{
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    GRBVar vars[n][n][n];
    int i, j, v;

    // Create 3-D array of model variables

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        for (v = 0; v < n; v++) {
          string s = "G_" + itos(i) + "_" + itos(j) + "_" + itos(v);
          vars[i][j][v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, s);
        }
      }
    }

    // Add constraints

    // Each cell must take one value

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        GRBLinExpr expr = 0;
        for (v = 0; v < n; v++)
          expr += vars[i][j][v];
        string s = "V_" + itos(i) + "_" + itos(j);
        model.addConstr(expr, GRB_EQUAL, 1.0, s);
      }
    }

    // Each value appears once per row

    for (i = 0; i < n; i++) {
      for (v = 0; v < n; v++) {
        GRBLinExpr expr = 0;
        for (j = 0; j < n; j++)
          expr += vars[i][j][v];
        string s = "R_" + itos(i) + "_" + itos(v);
        model.addConstr(expr == 1.0, s);
      }
    }

    // Each value appears once per column

    for (j = 0; j < n; j++) {
      for (v = 0; v < n; v++) {
        GRBLinExpr expr = 0;
        for (i = 0; i < n; i++)
          expr += vars[i][j][v];
        string s = "C_" + itos(j) + "_" + itos(v);
        model.addConstr(expr == 1.0, s);
      }
    }

    // Each value appears once per sub-grid

    for (v = 0; v < n; v++) {
      for (int i0 = 0; i0 < sd; i0++) {
        for (int j0 = 0; j0 < sd; j0++) {
          GRBLinExpr expr = 0;
          for (int i1 = 0; i1 < sd; i1++) {
            for (int j1 = 0; j1 < sd; j1++) {
              expr += vars[i0*sd+i1][j0*sd+j1][v];
            }
          }

          string s = "Sub_" + itos(v) + "_" + itos(i0) + "_" + itos(j0);
          model.addConstr(expr == 1.0, s);
        }
      }
    }

    // Fix variables associated with pre-specified cells

    char input[10];
    for (i = 0; i < n; i++) {
      cin >> input;
      for (j = 0; j < n; j++) {
        int val = (int) input[j] - 48 - 1; // 0-based

        if (val >= 0)
          vars[i][j][val].set(GRB_DoubleAttr_LB, 1.0);
      }
    }

    // Optimize model

    model.optimize();

    // Write model to file

    model.write("sudoku.lp");

    cout << endl;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        for (v = 0; v < n; v++) {
          if (vars[i][j][v].get(GRB_DoubleAttr_X) > 0.5)
            cout << v+1;
        }
      }
      cout << endl;
    }
    cout << endl;
  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Error during optimization" << endl;
  }

  return 0;
}

/* Copyright 2024, Gurobi Optimization, LLC */

/* This example formulates and solves the following simple MIP model:

     maximize    x +   y + 2 z
     subject to  x + 2 y + 3 z <= 4
                 x +   y       >= 1
                 x, y, z binary
*/

using namespace std;

int
main(int   argc,
     char *argv[])
{
  try {

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
    GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
    GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

    // Set objective: maximize x + y + 2 z
    model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

    // Add constraint: x + 2 y + 3 z <= 4
    model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

    // Add constraint: x + y >= 1
    model.addConstr(x + y >= 1, "c1");

    // Optimize model
    model.optimize();

    cout << x.get(GRB_StringAttr_VarName) << " "
         << x.get(GRB_DoubleAttr_X) << endl;
    cout << y.get(GRB_StringAttr_VarName) << " "
         << y.get(GRB_DoubleAttr_X) << endl;
    cout << z.get(GRB_StringAttr_VarName) << " "
         << z.get(GRB_DoubleAttr_X) << endl;

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}

/* Copyright 2024, Gurobi Optimization, LLC */

/* This example reads a MIP model from a file, solves it and
   prints the objective values from all feasible solutions
   generated while solving the MIP. Then it creates the fixed
   model and solves that model. */

//#include <cmath>
using namespace std;

int
main(int   argc,
     char *argv[])
{
  if (argc < 2) {
    cout << "Usage: mip2_c++ filename" << endl;
    return 1;
  }

  GRBEnv *env = 0;
  GRBVar *fvars = 0;
  try {
    env = new GRBEnv();
    GRBModel model = GRBModel(*env, argv[1]);

    if (model.get(GRB_IntAttr_IsMIP) == 0) {
      throw GRBException("Model is not a MIP");
    }

    model.optimize();

    int optimstatus = model.get(GRB_IntAttr_Status);

    cout << "Optimization complete" << endl;
    double objval = 0;
    if (optimstatus == GRB_OPTIMAL) {
      objval = model.get(GRB_DoubleAttr_ObjVal);
      cout << "Optimal objective: " << objval << endl;
    } else if (optimstatus == GRB_INF_OR_UNBD) {
      cout << "Model is infeasible or unbounded" << endl;
      return 0;
    } else if (optimstatus == GRB_INFEASIBLE) {
      cout << "Model is infeasible" << endl;
      return 0;
    } else if (optimstatus == GRB_UNBOUNDED) {
      cout << "Model is unbounded" << endl;
      return 0;
    } else {
      cout << "Optimization was stopped with status = "
           << optimstatus << endl;
      return 0;
    }

    /* Iterate over the solutions and compute the objectives */

    cout << endl;
    for ( int k = 0; k < model.get(GRB_IntAttr_SolCount); ++k ) {
      model.set(GRB_IntParam_SolutionNumber, k);
      double objn = model.get(GRB_DoubleAttr_PoolObjVal);

      cout << "Solution " << k << " has objective: " << objn << endl;
    }
    cout << endl;

    /* Create a fixed model, turn off presolve and solve */

    GRBModel fixed = model.fixedModel();

    fixed.set(GRB_IntParam_Presolve, 0);

    fixed.optimize();

    int foptimstatus = fixed.get(GRB_IntAttr_Status);

    if (foptimstatus != GRB_OPTIMAL) {
      cerr << "Error: fixed model isn't optimal" << endl;
      return 0;
    }

    double fobjval = fixed.get(GRB_DoubleAttr_ObjVal);

    if (fabs(fobjval - objval) > 1.0e-6 * (1.0 + fabs(objval))) {
      cerr << "Error: objective values are different" << endl;
      return 0;
    }

    int numvars = model.get(GRB_IntAttr_NumVars);

    /* Print values of nonzero variables */
    fvars = fixed.getVars();
    for (int j = 0; j < numvars; j++) {
      GRBVar v = fvars[j];
      if (v.get(GRB_DoubleAttr_X) != 0.0) {
        cout << v.get(GRB_StringAttr_VarName) << " "
             << v.get(GRB_DoubleAttr_X) << endl;
      }
    }

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Error during optimization" << endl;
  }

  delete[] fvars;
  delete env;
  return 0;
}