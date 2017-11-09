#include "ruby.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "vtc.h"
#include "vtclocal.h"

static double mmin = 0.0;
static double xmin = 0.0, xmax = 0.0;
static int autoscale = 1;
static int node_div_crit = 8;
static double eps = 0.02;
static double theta = 0.75;
static int ncrit = 100;
static int grape_flag = 0;
static double pprad = 0.2;
static int full_dof_flag = TRUE;
static int me_order = 1;
static int negativemass = FALSE;

typedef struct {
  int n;
  double *airdrag;
  double *admas;
  double *length;
  double *hook;
  long (*edge)[2];
} Springraph;

VALUE cNbody;



static char ptype(VALUE obj)
{
  switch (TYPE(obj)) {
  case T_FLOAT:
    printf ("float\n");
    break;
  case T_FIXNUM:
    printf ("fixnum\n");
    break;
  case T_STRING:
    printf ("string\n");
    break;
  case T_ARRAY:
    printf ("array\n");
    break;
  case T_HASH:
    printf ("hash\n");
    break;
  case T_STRUCT:
    printf ("struct\n");
    break;
  case T_BIGNUM:
    printf ("bignum\n");
    break;
  case T_TRUE:
    printf ("true\n");
    break;
  case T_FALSE:
    printf ("false\n");
    break;
  case T_DATA:
    printf ("data\n");
    break;
  case  T_MATCH:
    printf ("match\n");
    break;
  case  T_SYMBOL:
    printf ("symbol\n");
    break;
  default:
    /* 例外を発生させる */
    rb_raise(rb_eTypeError, "unknown.");
    break;
  }
}

create_initcond(Nbodyinfo *nb, int n)
{
    int i, k;
    double m = 1.0/n;

    nb->n = n;
    nb->m = (double *)malloc(sizeof(double)*n);
    if (NULL == nb->m) {
	perror("create_initcond");
	exit(1);
    }
    nb->x = (double (*)[3])malloc(sizeof(double)*3*n);
    nb->v = (double (*)[3])malloc(sizeof(double)*3*n);
    if (NULL == nb->x || NULL == nb->v) {
	perror("create_initcond");
	exit(1);
    }
    for (i = 0; i < n; i++) {
	nb->m[i] = m;
    }
    for (i = 0; i < n; i++) {
	nb->v[i][0] = 0.0;
	nb->v[i][1] = 0.0;
	nb->v[i][2] = 0.0;
	nb->a[i][0] = 0.0;
	nb->a[i][1] = 0.0;
	nb->a[i][2] = 0.0;
    }
    i = 0;
    while (i < n) {
	double x, y, z;
	x = 2.0*drand48()-1.0;
	y = 2.0*drand48()-1.0;
	z = 2.0*drand48()-1.0;
	if (x*x+y*y+z*z<1.0) {
	    nb->x[i][0] = x;
	    nb->x[i][1] = y;
	    nb->x[i][2] = z;
	    i++;
	}
    }
}

static void
init_nbodyinfo(Nbodyinfo *nb)
{
    int k;

    nb->m = NULL;
    nb->x = NULL;
    nb->v = NULL;
    nb->a = NULL;
    nb->p = NULL;
}


void free_forceinfo(Forceinfo *fi)
{
  if (fi)
    {
    free(fi);
  }
}

void free_nbodyinfo(Nbodyinfo *nb)
{
  if (nb->m)
    {
    free(nb->m);
  }
  if (nb->x)
    {
    free(nb->x);
  }
  if (nb->v)
    {
    free(nb->v);
  }
  if (nb->a)
    {
    free(nb->a);
  }
  if (nb)
    {
    free(nb);
  }
#ifdef WITHGRAPE
  grape_close();
#endif
}

static void
print_nb(Nbodyinfo *nb)
{
    int i;
    int n = nb->n;
    double *nbx0, *nbx1, *nbx2; 
    double *nbv0, *nbv1, *nbv2; 

    nbx0 = (double *)&(nb->x[0][0]); 
    nbx1 = (double *)&(nb->x[0][1]); 
    nbx2 = (double *)&(nb->x[0][2]); 
    nbv0 = (double *)&(nb->v[0][0]); 
    nbv1 = (double *)&(nb->v[0][1]); 
    nbv2 = (double *)&(nb->v[0][2]); 
    for (i = 0; i < n; i++)
      { 
	fprintf(stdout, "i:%d: ", i); 
	fprintf(stdout, "position %f, ", nbx0[i*3]); 
	fprintf(stdout, "%f, ", nbx1[i*3]); 
	fprintf(stdout, "%f\n", nbx2[i*3]); 
	fprintf(stdout, "velocity %f, ", nbv0[i*3]); 
	fprintf(stdout, "%f, ", nbv1[i*3]); 
	fprintf(stdout, "%f\n", nbv2[i*3]); 
    } 
}

static VALUE nb_grape_close(VALUE self)
{
  grape_close();
}



static VALUE nb_get_mass(VALUE self)
{
  VALUE ret, np;
  Nbodyinfo *nb;
  int n, i;

  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  n = nb->n;
  ret = rb_ary_new();
  for (i = 0; i < n; i++)  
    {
      rb_ary_push(ret,  rb_float_new(nb->m[i]));
    }
  return ret;  
}

static VALUE nb_get_position(VALUE self)
{
  VALUE np, ret, item;
  Nbodyinfo *nb;
  int n, i, j;

  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  n = nb->n;
  ret = rb_ary_new();
  for (i = 0; i < n; i++)  
    {
      item = rb_ary_new();
      for (j = 0; j < 3; j++)  
 	{  
  	  rb_ary_push(item,  rb_float_new(nb->x[i][j]));
  	}
      rb_ary_push(ret, item);
    }
  return ret;
}


static VALUE nb_get_velocity(VALUE self)
{
  VALUE np, ret, item;
  Nbodyinfo *nb;
  int n, i, j;

  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  n = nb->n;
  ret = rb_ary_new();
  for (i = 0; i < n; i++)  
    {
      item = rb_ary_new();
      for (j = 0; j < 3; j++)  
 	{  
  	  rb_ary_push(item,  rb_float_new(nb->v[i][j]));
  	}
      rb_ary_push(ret, item);
    }
  return ret;
}

static VALUE nb_get_acc(VALUE self)
{
  VALUE np, pos, pitem;
  Nbodyinfo *nb;
  int n, i, j;

  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  n = nb->n;
  pos = rb_ary_new();
  for (i = 0; i < n; i++)  
    {
      pitem = rb_ary_new();
      for (j = 0; j < 3; j++)  
 	{  
  	  rb_ary_push(pitem,  rb_float_new(nb->a[i][j]));
  	}
      rb_ary_push(pos, pitem);
    }
  return pos;
}

static VALUE nb_read_mass(VALUE self, VALUE input)
{
  VALUE np, item;
  Nbodyinfo *nb;
  int i;
  long n;

  n = RARRAY_LEN(input); 
  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  for (i = 0; i < n; i++)
    {  
      nb->m[i] = rb_float_value(rb_ary_entry(input, i));
    }
  return self;
}

static VALUE nb_read_position(VALUE self, VALUE input)
{
  VALUE np, item;
  Nbodyinfo *nb;
  int i, j;
  long n;

  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  n = RARRAY_LEN(input); 
  nb->n = n;
  for (i = 0; i < n; i++)
    {  
      item = rb_ary_entry(input, i);	  
      for (j = 0; j < 3; j++)
	{
	  nb->x[i][j] = rb_float_value(rb_ary_entry(item, j));
      }
    }
  return self;
}


static VALUE nb_read_velocity(VALUE self, VALUE input)
{
  VALUE np, item;
  Nbodyinfo *nb;
  int i, j;
  long n;

  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  n = RARRAY_LEN(input); 
  for (i = 0; i < n; i++)
    {  
      item = rb_ary_entry(input, i);	  
      for (j = 0; j < 3; j++)
	{
	  nb->v[i][j] = rb_float_value(rb_ary_entry(item, j));
      }
    }
  return self;
}

static VALUE nb_read_acc(VALUE self, VALUE input)
{
  VALUE np, item;
  Nbodyinfo *nb;
  int i, j;
  long n;

  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  n = RARRAY_LEN(input); 
  for (i = 0; i < n; i++)
    {  
      item = rb_ary_entry(input, i);	  
      for (j = 0; j < 3; j++)
	{
	  nb->a[i][j] = rb_float_value(rb_ary_entry(item, j));
	}
    }
  return self;
}

static void push_velocity(double dt, Nbodyinfo *nb)
{
  int i, j;
  int n = nb->n;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  nb->v[i][j] += dt*(nb->a[i][j]);
	}
    }
}

static void push_velocity_dg(VALUE self, double dt, Nbodyinfo *nb)
{
  int i, j;
  int n = nb->n;
  Springraph *sp;
  VALUE sip;

  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);

  for (i = 0; i < n; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  nb->v[i][j] += dt*(nb->a[i][j] - sp->airdrag[i]*(nb->v[i][j]));
	}
    }
}

static void push_position(double dt, Nbodyinfo *nb)
{
  int i, j;
  int n = nb->n;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  nb->x[i][j] += dt*(nb->v[i][j]);
	}
    }
}

static VALUE nb_frog(VALUE self, VALUE dt)
{
  Nbodyinfo *nb;
  VALUE np;
  float slice;

  slice = rb_float_value(dt);
  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  push_velocity(0.5*slice, nb);
  push_position(slice, nb);
  push_velocity(0.5*slice, nb);
  return self;
}

static VALUE nb_dg_frog(VALUE self, VALUE dt)
{
  Nbodyinfo *nb;
  VALUE np;
  float slice;

  slice = rb_float_value(dt);
  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  push_velocity_dg(self, 0.5*slice, nb);
  push_position(slice, nb);
  push_velocity_dg(self, 0.5*slice, nb);
  return self;
}

static void free_edgeinfo(Springraph *sp)
{
  if (sp->admas)
    {
    free(sp->admas);
  }
  if (sp->length)
    {
    free(sp->length);
  }
  if (sp->edge)
    {
    free(sp->edge);
  }
  if (sp->airdrag)
    {
    free(sp->airdrag);
  }
  if (sp)
    {
    free(sp);
  }

}

static VALUE nb_load_edge(VALUE self, VALUE edge)
{
  Springraph *sp;
  Nbodyinfo *nb;
  VALUE sip, item, np;
  long n, i, j;
  int vn;

  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  vn=nb->n;

  n = RARRAY_LEN(edge);

  sp = (Springraph *)malloc(sizeof(Springraph));
  sp->edge = (long (*)[2])malloc(sizeof(long)*2*n);
  sp->hook = (double *)malloc(sizeof(double)*n);
  sp->admas = (double *)malloc(sizeof(double)*vn);
  sp->airdrag = (double *)malloc(sizeof(double)*vn);

  sp->n = n;
  for (i = 0; i < n; i++)
    {
      item = rb_ary_entry(edge, i);
      for (j = 0; j < 2; j++)
	{
	  sp->edge[i][j] = FIX2INT(rb_ary_entry(item, j));
      }
    }
  sip = Data_Wrap_Struct(cNbody, 0, free_edgeinfo, sp);
  rb_iv_set(self, "@si", sip);
  return self;
}

static VALUE nb_get_hookparams(VALUE self)
{
  VALUE ret, np, sip;
  Springraph *sp;
  int n, i;

  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);
  n = sp->n;
  ret = rb_ary_new();
  for (i = 0; i < n; i++)  
    {
      rb_ary_push(ret,  rb_float_new(sp->hook[i]));
    }
  return ret;
}

static VALUE nb_load_hookparams(VALUE self, VALUE input)
{
  Springraph *sp;
  VALUE sip;
  int i, j, n;

  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);

  n = sp->n;
  for (i = 0; i < n; i++)
    {
      sp->hook[i] = rb_float_value(rb_ary_entry(input, i));
    }
  return self;
}

static VALUE nb_get_additional_mass(VALUE self)
{
  Springraph *sp;
  Nbodyinfo *nb;
  VALUE sip, np, ret;
  int i, j, n;

  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);
  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);

  n = nb->n;
  ret = rb_ary_new();
  for (i = 0; i < n; i++)
    {
      rb_ary_push(ret,  rb_float_new(sp->admas[i]));
    }
  return ret;
}


static VALUE nb_load_additional_mass(VALUE self, VALUE input)
{
  Springraph *sp;
  Nbodyinfo *nb;
  VALUE sip, np;
  int i, j, n;

  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);
  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);

  n = nb->n;
  for (i = 0; i < n; i++)
    {
      sp->admas[i] = rb_float_value(rb_ary_entry(input, i));
    }
  return self;
}


static VALUE nb_get_airdrag(VALUE self)
{
  Springraph *sp;
  Nbodyinfo *nb;
  VALUE sip, np, ret;
  int i, j, n;

  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);
  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);

  n = nb->n;
  ret = rb_ary_new();
  for (i = 0; i < n; i++)
    {
      rb_ary_push(ret,  rb_float_new(sp->airdrag[i]));
    }
  return ret;
}


static VALUE nb_load_airdrag(VALUE self, VALUE input)
{
  Springraph *sp;
  Nbodyinfo *nb;
  VALUE sip, np;
  int i, j, n;

  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);
  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);

  n = nb->n;
  for (i = 0; i < n; i++)
    {
      sp->airdrag[i] = rb_float_value(rb_ary_entry(input, i));
    }
  return self;
}

static VALUE nb_get_spring_length(VALUE self)
{
  VALUE ret, np, sip;
  Springraph *sp;
  int n, i;

  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);
  n = sp->n;
  ret = rb_ary_new();
  for (i = 0; i < n; i++)  
    {
      rb_ary_push(ret,  rb_float_new(sp->length[i]));
    }
  return ret;
}

static VALUE nb_load_spring_length(VALUE self, VALUE input)
{
  Springraph *sp;
  VALUE sip;
  int i, j, n;

  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);

  n = sp->n;
  sp->length = (double *)malloc(sizeof(double)*n);
  for (i = 0; i < n; i++)
    {
      sp->length[i] = rb_float_value(rb_ary_entry(input, i));
    }
  return self;
}

static VALUE nb_edge_list(VALUE self)
{
  Springraph *sp;
  VALUE sip, np, ret, item;
  int i, j, n, p, q;
  double h, l, d;
  
  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);
  n = sp->n;
  ret = rb_ary_new();
  for (i = 0; i < n; i++)  
    {
      item = rb_ary_new();
      for (j = 0; j < 2; j++)  
 	{  
  	  rb_ary_push(item,  INT2FIX(sp->edge[i][j]));
  	}
      rb_ary_push(ret, item);
    }
  return ret;  
}

static VALUE nb_add_magnetic(VALUE self, VALUE magn, VALUE direct)
{
  Springraph *sp;
  Nbodyinfo *nb;
  VALUE sip, np;
  double mg, d;
  double dir[3], ev[3];
  int i, j, p, q, n;

  mg = rb_float_value(magn);
  for (i = 0; i < 2; i++)
    {
      dir[i] = rb_float_value(rb_ary_entry(direct, i));
    }
  
  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);
  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);

  n = sp->n;
  for (i = 0; i < n; i++)
    {
      p = sp->edge[i][1];
      q = sp->edge[i][0];

      d = 0.0;
      for (j = 0; j < 3; j++)
	{
	  ev[j] = nb->x[p][j] - nb->x[q][j];
	  d = d + (ev[j])*(ev[j]);
	}
      d = sqrt(d);

      if (d > 0) {
	for (j = 0; j < 3; j++)
	  {
	    nb->a[p][j] = nb->a[p][j] + (mg * (dir[j] - (ev[j] / d)))/sp->admas[p];
	    nb->a[q][j] = nb->a[q][j] - (mg * (dir[j] - (ev[j] / d)))/sp->admas[q];
	  }
      }
    }
  return self;
}

static VALUE nb_add_spring(VALUE self)
{
  Springraph *sp;
  Nbodyinfo *nb;
  VALUE sip, np;
  int i, j, n, p, q, vn;
  double h, l, m, d, a0, a1, ha;
  
  sip = rb_iv_get(self, "@si");
  Data_Get_Struct(sip, Springraph, sp);
  np = rb_iv_get(self, "@ni");
  Data_Get_Struct(np, Nbodyinfo, nb);
  
  n = sp->n;
  vn = nb->n;
  for (i = 0; i < n; i++)
    {
      p = sp->edge[i][0];
      q = sp->edge[i][1];
      h = sp->hook[i];
      l = sp->length[i];

      d = 0.0;
      for (j = 0; j < 3; j++) {
	d = d + ((nb->x[p][j] - nb->x[q][j]) * (nb->x[p][j] - nb->x[q][j]));}
      d = sqrt(d);
      
      if (d > 0) {
	for (j = 0; j < 3; j++)
	  {
	    ha = (h * ((nb->x[p][j] - nb->x[q][j]) * (d - l)) / d);
	    nb->a[p][j] = nb->a[p][j] - ha/sp->admas[p];
	    nb->a[q][j] = nb->a[q][j] + ha/sp->admas[q];
	  }
      }
    }
  return self;
}

static VALUE nb_get_bhvector_param(VALUE self)
{
  Forceinfo *fi;
  VALUE fp;

  fp = rb_iv_get(self, "@fi");
  Data_Get_Struct(fp, Forceinfo, fi);

  return INT2FIX(fi->ncrit);
}

static VALUE nb_set_bhvector_param(VALUE self, VALUE input)
{
  Forceinfo *fi;
  VALUE fp;
  
  fp = rb_iv_get(self, "@fi");
  Data_Get_Struct(fp, Forceinfo, fi);

  fi->ncrit = FIX2INT(input);
  return self;
}

  

static VALUE nb_init(VALUE self, VALUE num)
/* static VALUE nb_init(VALUE self) */
{
  Forceinfo *fi;
  Nbodyinfo *ni;
  Springraph *sp;
  VALUE fp, np, ah, ary;
  Nbodyinfo *t;
  int n = 3;

  n = NUM2INT(num);
  fi = (Forceinfo *)malloc(sizeof(Forceinfo));
  ni = (Nbodyinfo *)malloc(sizeof(Nbodyinfo));

  init_nbodyinfo(ni);
  (*ni).a = (double (*)[3])malloc(sizeof(double)*3*n);
  (*ni).p = (double *)malloc(sizeof(double)*n);
  create_initcond(ni, n);
  
  vtc_get_default_tree_params(fi);
/* #ifdef WITHGRAPE */
  fi->calculator = GRAPE;
/* #endif */

  /* fprintf(stdout, "before %d\n", &nb); */
  fp = Data_Wrap_Struct(cNbody, 0, free_forceinfo, fi);
  np = Data_Wrap_Struct(cNbody, 0, free_nbodyinfo, ni);
  rb_iv_set(self, "@fi", fp);
  rb_iv_set(self, "@ni", np);

/*   ah = rb_iv_get(self, "@fi"); */
/*   Data_Get_Struct(ah, Forceinfo, t);   */

  return self;
}


static VALUE nbody_get_force_tree(VALUE self)
{
  Forceinfo *fi;
  Nbodyinfo *nb;
  Springraph *sb;
  VALUE np, fp, sp;
  int n, i, j;

  fp = rb_iv_get(self, "@fi");
  np = rb_iv_get(self, "@ni");
  sp = rb_iv_get(self, "@si");

  Data_Get_Struct(fp, Forceinfo, fi);
  Data_Get_Struct(np, Nbodyinfo, nb);
  Data_Get_Struct(sp, Springraph, sb);

  n = nb->n;

  vtc_get_force_tree(fi, nb);
  return self;
}

static VALUE nbody_get_coulomb_force(VALUE self)
{
  Forceinfo *fi;
  Nbodyinfo *nb;
  Springraph *sb;
  VALUE np, fp, sp;
  int n, i, j;

  fp = rb_iv_get(self, "@fi");
  np = rb_iv_get(self, "@ni");
  sp = rb_iv_get(self, "@si");

  Data_Get_Struct(fp, Forceinfo, fi);
  Data_Get_Struct(np, Nbodyinfo, nb);
  Data_Get_Struct(sp, Springraph, sb);

  n = nb->n;

  vtc_get_force_tree(fi, nb);

  for (i=0; i < n ; i++)
    {
      for (j=0; j < 3; j++)
  	{ 
  	  nb->a[i][j] = - nb->a[i][j] * nb->m[i] / sb->admas[i];
  	}
    }
  
  return self;
}


void Init_Nbody() {
  int n;

  cNbody = rb_define_class("Nbody", rb_cObject);
/*   rb_define_method(cNbody, "initialize", nb_init, 0); */
  rb_define_method(cNbody, "initialize", nb_init, 1);
  rb_define_method(cNbody, "force", nbody_get_force_tree, 0);
  rb_define_method(cNbody, "coulomb_force", nbody_get_coulomb_force, 0);
  rb_define_method(cNbody, "chg", nb_get_mass, 0);
  rb_define_method(cNbody, "pos", nb_get_position, 0);
  rb_define_method(cNbody, "vel", nb_get_velocity, 0);
  rb_define_method(cNbody, "acc", nb_get_acc, 0);
  rb_define_method(cNbody, "chg=", nb_read_mass, 1);
  rb_define_method(cNbody, "pos=", nb_read_position, 1);
  rb_define_method(cNbody, "vel=", nb_read_velocity, 1);
  rb_define_method(cNbody, "acc=", nb_read_acc, 1);
  rb_define_method(cNbody, "frog", nb_dg_frog, 1);
  rb_define_method(cNbody, "pure_frog", nb_frog, 1);
  rb_define_method(cNbody, "edge=", nb_load_edge, 1);
  rb_define_method(cNbody, "edge", nb_edge_list, 0);
  rb_define_method(cNbody, "bane", nb_add_spring, 0);
  rb_define_method(cNbody, "mgn", nb_add_magnetic, 2);
  rb_define_method(cNbody, "springlength=", nb_load_spring_length, 1);
  rb_define_method(cNbody, "springlength", nb_get_spring_length, 0);
  rb_define_method(cNbody, "mas=", nb_load_additional_mass, 1);
  rb_define_method(cNbody, "mas", nb_get_additional_mass, 0);
  rb_define_method(cNbody, "hook=", nb_load_hookparams, 1);
  rb_define_method(cNbody, "hook", nb_get_hookparams, 0);
  rb_define_method(cNbody, "tree_vector_param", nb_get_bhvector_param, 0);
  rb_define_method(cNbody, "tree_vector_param=", nb_set_bhvector_param, 1);
  rb_define_method(cNbody, "airdrag=", nb_load_airdrag, 1);
  rb_define_method(cNbody, "airdrag", nb_get_airdrag, 0);
#ifdef WITHGRAPE
  rb_define_singleton_method(cNbody, "grape_close", nb_grape_close, 0);
#endif
}
