#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <errno.h>
#include <getopt.h>
#include <stdint.h>

#include <vector>
#include <string>

double _drand() {
  return (double)(rand())/(RAND_MAX + 1.0);
}

typedef struct sp_type {
  int32_t m_nclause, m_nvar;

  // clauses with each variable, packed.
  // Note: variables start at 1 and use negative values to indicate a negation.
  // example:
  //
  //    ... 1 -5 4 -1 5 3 4 ...
  //
  // could represent two clauses:  (1 + -5 + 4) . (-1 + 5 + 3 + 4)
  // m_clause_info has the information on where each clause starts and it's length
  //
  std::vector<int32_t> m_clause;

  // interleaved start of clause, size.
  // example:
  //
  //    ... 217 3 ...
  //
  // means clause starts at m_clause[217] and is 3 long
  //
  std::vector<int32_t> m_clause_info;

  // where each variable appears in the clause,
  // interleaved with index
  //
  std::vector<int32_t> m_var;

  // interleaved start and size of variable data
  // example:
  //
  //    ... 31 7 ...
  //
  // means variable information starts at m_var[31] and is 7 long
  //
  std::vector<int32_t> m_var_info;


  int32_t m_t;
  std::vector<double> m_c2v[2];
  std::vector<double> m_v2c_u[2];
  std::vector<double> m_v2c_s[2];
  std::vector<double> m_v2c_x[2];

  sp_type() {
    m_nclause = 0;
    m_nvar = 0;
  }

  void init_sp() {
    int32_t i, n;
    m_c2v[0].clear();
    m_c2v[1].clear();

    m_v2c_u[0].clear();
    m_v2c_u[1].clear();

    m_v2c_s[0].clear();
    m_v2c_s[1].clear();

    m_v2c_x[0].clear();
    m_v2c_x[1].clear();

    m_t = 0;
    n = m_clause.size()/2;
    for (i=0; i<n; i++) {
      m_c2v[0].push_back( _drand() );
      m_c2v[1].push_back( 0.0 );
    }

    n = m_var.size()/2;
    for (i=0; i<n; i++) {
      m_v2c_u[0].push_back( _drand() );
      m_v2c_s[0].push_back( _drand() );
      m_v2c_x[0].push_back( _drand() );

      m_v2c_u[1].push_back( 0.0 );
      m_v2c_s[1].push_back( 0.0 );
      m_v2c_x[1].push_back( 0.0 );
    }

  }

  //TODO:
  // C(i) = { a \in C : i \in V(a) }
  //   - clauses set of variable i
  //
  // C_a^s(i) = { b \in C(i) \ a : J_{a,i} = J_{b,i}
  //   - clause set that have same sign as var i
  // C_a^u(i) = { b \in C(i) \ a : J_{a,i} \ne J_{b,i}
  //   - clause set that have diff. sign as var i
  //
  // implement lookup
  //

  void sp_tick_v2c(void) {
    int32_t t_cur, t_nxt;
    int32_t i, j, k, n;

    int32_t clause, var;

    double prod_sat, prod_unsat;

    t_cur = m_t;
    t_nxt = 1-m_t;

    for (clause=0; clause < m_nclause; clause++) {

      c_s = m_clause_info[2*clause];
      c_n = m_clause_info[2*clause+1];

      for (c_idx=0; c_idx < c_n; c_idx++) {
        varname = m_clause[c_s + c_idx];
        if (varname<0) { varname = -varname; }

      }
    }



  }

  void sp_tick_var(void) {
  }

  void print_dimacs(FILE *fp) {
    int32_t i, j, p, n;

    fprintf(fp, "p %i %i\n", (int)m_nvar, (int)m_nclause);
    for (i=0; i<m_clause_info.size(); i+=2) {
      p = m_clause_info[i];
      n = m_clause_info[i+1];
      for (j=0; j<n; j++) {
        printf("%i", (int)(m_clause[p+j]));
      }
      printf("0\n");
    }
  }

  void print_sp(void) {
    int32_t i, j, p, n;
    int32_t _v, _idx;

    printf("---\n\n");
    printf("m_nclause: %i, m_nvar: %i\n", (int)m_nclause, (int)m_nvar);
    printf("\n---\n");

    for (i=0; i<m_clause_info.size(); i+=2) {
      p = m_clause_info[i];
      n = m_clause_info[i+1];
      printf("c%i [%i+%i]{%i}:", (int)(i/2), (int)p, (int)n, (int)m_t);
      for (j=0; j<n; j++) {
        printf(" %i {%f}", (int)(m_clause[p+j]), (float)m_c2v[m_t][p+j]);
      }
      printf("\n");
    }

    printf("\n---\n");

    for (i=0; i<m_var_info.size(); i+=2) {

      _v = i/2;

      p = m_var_info[i];
      n = m_var_info[i+1];
      printf("v%i [%i+%i]{%i}:", (int)((i/2)+1), (int)p, (int)n, (int)m_t);
      for (j=0; j<n; j+=2) {

        _idx = j/2;

        printf(" c%i[%i]{%f,%f,%f}", m_var[p + j], m_var[p+j+1],
            (float)m_v2c_u[m_t][_v+_idx],
            (float)m_v2c_s[m_t][_v+_idx],
            (float)m_v2c_x[m_t][_v+_idx] );
      }
      printf("\n");
    }

    printf("---\n");

  }

  void print(void) {
    int32_t i, j, p, n;

    printf("---\n\n");
    printf("m_nclause: %i, m_nvar: %i\n", (int)m_nclause, (int)m_nvar);
    printf("\n---\n");

    for (i=0; i<m_clause_info.size(); i+=2) {
      p = m_clause_info[i];
      n = m_clause_info[i+1];
      printf("c%i [%i+%i]:", (int)(i/2), (int)p, (int)n);
      for (j=0; j<n; j++) {
        printf(" %i", (int)(m_clause[p+j]));
      }
      printf("\n");
    }

    printf("\n\n---\n");

    printf("m_clause[%i]:\n", (int)m_clause.size());
    for (i=0; i<m_clause.size(); i++) {
      if ((i%20)==0) { printf("\n"); }
      printf(" %i", m_clause[i]);
    }
    printf("\n\n");
    printf("---\n");

    for (i=0; i<m_var_info.size(); i+=2) {
      p = m_var_info[i];
      n = m_var_info[i+1];
      printf("v%i [%i+%i]:", (int)((i/2)+1), (int)p, (int)n);
      for (j=0; j<n; j+=2) {
        printf(" c%i[%i]", m_var[p + j], m_var[p+j+1]);
      }
      printf("\n");
    }

    printf("---\n");

    printf("m_var[%i]:\n", (int)m_var.size());
    for (i=0; i<m_var.size(); i++) {
      if ((i%20)==0) { printf("\n"); }
      printf(" %2i", m_var[i]);
    }
    printf("\n\n");

    printf("---\n");

  }
} sp_t;

struct option g_long_opt [] = {
  {0, 0, 0, 0},
};

void show_help(FILE *fp) {
  fprintf(fp, "\n\nusage:\n");
  fprintf(fp, "  survey-prop [-h]\n");
  fprintf(fp, "\n");
  fprintf(fp, "  [-h]     show help (this screen)\n");
}

int parse_int_line(std::string &str, std::vector<int32_t> &v) {
  size_t n;
  int32_t ival;
  int p;
  std::string tok;

  n = str.size();
  for (p=0; p<n; p++) {
    if ((str[p] == ' ') || 
        (str[p] == '\n') ||
        (str[p] == '\r')) {
      if (tok.size()>0) {
        ival = (int32_t)strtol(tok.c_str(), NULL, 10);
        v.push_back(ival);
      }
      tok.clear();
      continue;
    }
    tok += str[p];
  }

  if (tok.size()>0) {
    ival = (int32_t)strtol(tok.c_str(), NULL, 10);
    v.push_back(ival);
  }


  return 0;
}

int read_dimacs(FILE *fp, sp_t &sp) {
  int ch, i;
  std::string buf;
  std::vector<int32_t> int_line;
  int32_t _start, _n;

  int32_t _v_n, _v_s, _v_idx;
  int32_t _c_n, _c_s, _c_idx;

  int32_t varname;

  std::vector<int32_t> v_nvar;

  while (!feof(fp)) {
    ch = fgetc(fp);
    if ((ch==EOF) || (ch=='\n')) {

      if ( (buf.size()>0) &&
           (((buf[0] == '-') ||
            ((buf[0] >= '0') && (buf[0] <= '9')))) ) {
        int_line.clear();
        parse_int_line(buf, int_line);

        _start = sp.m_clause.size();
        _n = 0;
        sp.m_clause_info.push_back(_start);
        for (i=0; i<int_line.size(); i++) {
          if (int_line[i] == 0) { break; }
          sp.m_clause.push_back(int_line[i]);

          if (abs(int_line[i]) >= sp.m_nvar) {
            sp.m_nvar = int_line[i];
          }
          _n++;
        }
        sp.m_clause_info.push_back(_n);
        sp.m_nclause++;

      }
      buf.clear();
      continue;
    }
    buf += ch;
  }

  if (buf.size()>0) {
    int_line.clear();
    parse_int_line(buf, int_line);

    _start = sp.m_clause.size();
    _n = 0;
    sp.m_clause_info.push_back(_start);
    for (i=0; i<int_line.size(); i++) {
      if (int_line[i] == 0) { break; }
      sp.m_clause.push_back(int_line[i]);

      if (abs(int_line[i]) >= sp.m_nvar) {
        sp.m_nvar = int_line[i];
      }
    }
    sp.m_clause_info.push_back(_n);
    sp.m_nclause++;
  }

  //---
  //---

  // parsing is done, do variable
  // hookup
  //

  // construct variable lookups
  //
  for (i=0; i<sp.m_nvar; i++) {
    v_nvar.push_back(0);
    sp.m_var_info.push_back(0);
    sp.m_var_info.push_back(0);
  }

  // get total variable count
  //
  _n=0;
  for (i=0; i<sp.m_clause.size(); i++) {
    varname = sp.m_clause[i];
    if (varname<0) { varname = -varname; }
    if (varname<=0) { return -1; }

    _v_s  = (varname-1)*2;
    sp.m_var_info[_v_s + 1]+=2;
    _n++;
  }

  // allocate and init m_var
  //
  for (i=0; i<_n; i++) {
    sp.m_var.push_back(0);
    sp.m_var.push_back(0);
  }

  // setup start pointers to m_var array
  //
  // m_var_info[ (varname-1) ]      = index in m_var
  // m_var_info[ (varname-1) + 1 ]  = length of m_var run
  //
  for (_n=0, i=0; i<sp.m_var_info.size(); i+=2) {
    sp.m_var_info[i] = _n;
    _n += sp.m_var_info[i+1];
  }

  // do another pass of clause information to setup
  // variable lookups
  //
  for (i=0; i<sp.m_clause_info.size(); i+=2) {
    _c_s = sp.m_clause_info[i];
    _c_n = sp.m_clause_info[i+1];

    for (_c_idx=0; _c_idx < _c_n; _c_idx++) {
      varname = sp.m_clause[_c_s + _c_idx];
      if (varname<0) { varname = -varname; }

      _v_s = sp.m_var_info[2*(varname-1)];
      _v_idx = v_nvar[varname-1];
      sp.m_var[_v_s + 2*_v_idx] = i/2;
      sp.m_var[_v_s + 2*_v_idx + 1] = _c_idx;
      v_nvar[varname-1]++;

    }

  }

  return 0;
}

int main(int argc, char **argv) {
  int ch, opt_idx, r;
  FILE *ifp = stdin;
  std::string ifn;

  sp_t sp;

  r=0;
  ifn = "-";

  while ((ch = getopt_long(argc, argv, "h", g_long_opt, &opt_idx)) >= 0) {
    switch(ch) {
      case 0:
        break;
      case 'h':
        show_help(stdout);
        exit(0);
        break;
      default:
        show_help(stderr);
        exit(-1);
        break;
    }
  }



  r = read_dimacs(ifp, sp);
  if (r<0) {
    perror(ifn.c_str());
  }

  sp.init_sp();

  sp.print();
  sp.print_sp();


}
