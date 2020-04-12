// License: CC0
//

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

typedef struct sp_edge_type {

  // 'c' claus
  // 'v' var
  //
  int32_t m_type;

  // clause name or variable name
  // m_id[0] - 'from'
  // m_id[1] - 'to'
  //
  int32_t m_id[2];

  // index, as opposed to id,
  // references the position in the m_E
  // array
  // m_idx[0] - 'from' index
  // m_idx[1] - 'to' index
  //
  int32_t m_idx[2];

  // Backpointer to connection edge entry.
  // Offset in array.
  //
  int32_t m_bp_offset;

  // -1 if variable appears un-negated
  //  1 if variable appears negated in clause
  //
  int32_t m_J;

  // clause to var - [0]
  // var to clause unsat - [0]
  // var to clause   sat - [1]
  // var to clause  star - [2]
  //
  double m_d[2][3];

} sp_edge_t;

typedef struct sp_type {
  double m_eps;

  int32_t m_nclause, m_nvar;

  std::vector< int32_t > m_raw_clause;
  std::vector< int32_t > m_raw_clause_info;

  std::vector< int32_t > m_v_idx, m_c_idx;
  std::vector< std::vector< sp_edge_t > > m_E;

  // time index {0,1}
  //
  int32_t m_t;

  // Variables are named starting at 1 but have
  //   structures that are 0 indexed.
  // Clauses are 0-indexed (start at 0)
  //

  sp_type() {
    m_nclause = 0;
    m_nvar = 0;
    m_t = 0;

    m_eps = 1.0 / (1024.0*1024.0*256.0);
  }

  void print_raw_clause() {
    int i;
    printf("nclause: %i, nvar: %i\n", m_nclause, m_nvar);
    for (i=0; i<m_raw_clause.size(); i++) { printf(" %i", m_raw_clause[i]); }
    printf("\n");
    for (i=0; i<m_raw_clause_info.size(); i++) { printf(" %i", m_raw_clause_info[i]); }
    printf("\n");
    printf("---\n");
  }

  int consistency(void) {
    int32_t i, j, k;
    int32_t c_idx, v_idx, p;
    int32_t _q, _s;

    if (m_E.size() != (m_nvar + m_nclause)) { return -100; }

    for (c_idx=0; c_idx < m_c_idx.size(); c_idx++) {
      p = m_c_idx[c_idx];

      for (i=0; i<m_E[p].size(); i++) {

        _q = m_E[p][i].m_idx[1];
        _s = m_E[p][i].m_bp_offset;

        if (m_E[p][i].m_id[0] != m_E[_q][_s].m_id[1]) { return -1; }
        if (m_E[p][i].m_id[1] != m_E[_q][_s].m_id[0]) { return -2; }

        if (m_E[p][i].m_idx[0] != m_E[_q][_s].m_idx[1]) { return -3; }
        if (m_E[p][i].m_idx[1] != m_E[_q][_s].m_idx[0]) { return -4; }

        if (m_E[p][i].m_J != m_E[_q][_s].m_J) { return -5; }

        if (m_E[p][i].m_type != 'c') { return -6; }

      }
    }

    for (v_idx=0; v_idx < m_v_idx.size(); v_idx++) {
      p = m_v_idx[v_idx];

      for (i=0; i<m_E[p].size(); i++) {
        if (m_E[p][i].m_type != 'v') { return -7; }

        _q = m_E[p][i].m_idx[1];
        _s = m_E[p][i].m_bp_offset;

        if (m_E[p][i].m_id[0] != m_E[_q][_s].m_id[1]) { return -8; }
        if (m_E[p][i].m_id[1] != m_E[_q][_s].m_id[0]) { return -9; }

        if (m_E[p][i].m_idx[0] != m_E[_q][_s].m_idx[1]) { return -10; }
        if (m_E[p][i].m_idx[1] != m_E[_q][_s].m_idx[0]) { return -11; }

        if (m_E[p][i].m_J != m_E[_q][_s].m_J) { return -12; }

      }
    }

    return 0;
  }


  void init_sp() {
    m_raw_clause.clear();
    m_raw_clause_info.clear();
    m_v_idx.clear();
    m_c_idx.clear();
    m_E.clear();
  }

  void print_E(void) {
    int32_t i, j, p;
    char _f, _t;

    for (i=0; i<m_E.size(); i++) {

      for (j=0; j<m_E[i].size(); j++) {
        printf("[%i][%i] %c ",
            i, j, (char)m_E[i][j].m_type);
        if      (m_E[i][j].m_type == (int)'c') { _f = 'c'; _t = 'v'; }
        else if (m_E[i][j].m_type == (int)'v') { _f = 'v'; _t = 'c'; }
        else                    { _f = 'x'; _t = 'x'; }

        printf("(%c%i -> %c%i) J%i {%f,%f,%f}{%f,%f,%f}\n",
            _f, m_E[i][j].m_id[0],
            _t, m_E[i][j].m_id[1],
            m_E[i][j].m_J,
            (float)m_E[i][j].m_d[0][0], (float)m_E[i][j].m_d[0][1], (float)m_E[i][j].m_d[0][2],
            (float)m_E[i][j].m_d[1][0], (float)m_E[i][j].m_d[1][1], (float)m_E[i][j].m_d[1][2]);
      }
    }
  }

  void print_dimacs(FILE *fp) {
    int32_t i, j, p, n;
    int32_t s;

    fprintf(fp, "p %i %i\n", (int)m_nvar, (int)m_nclause);
    for (i=0; i<m_c_idx.size(); i++) {
      p = m_c_idx[i];
      for (j=0; j<m_E[p].size(); j++) {
        if (j>0) { fprintf(fp, " "); }

        if (m_E[p][j].m_J < 0) {
          fprintf(fp, "%i", m_E[p][j].m_id[1]);
        }
        else {
          fprintf(fp, "-%i", m_E[p][j].m_id[1]);
        }
      }
      if (j>0) { fprintf(fp, " "); }
      fprintf(fp, "0\n");
    }

  }

  void print_dot(FILE *fp) {
    int32_t i, j, p, n;
    int32_t s;

    fprintf(fp, "graph G {\n");
    for (i=0; i<m_c_idx.size(); i++) {
      p = m_c_idx[i];


      for (j=0; j<m_E[p].size(); j++) {
        if (j==0) {
          fprintf(fp, "  c%i [shape=box];\n", m_E[p][j].m_id[0]);
        }

        fprintf(fp, "  c%i -- v%i", m_E[p][j].m_id[0], m_E[p][j].m_id[1]);
        if (m_E[p][j].m_J > 0) {
          fprintf(fp, " [style=dotted]");
        }
        fprintf(fp, ";\n");
      }
    }
    fprintf(fp, "}\n");
  }

  void print_dot_w(FILE *fp) {
    int32_t i, j, p, n;
    int32_t s;

    fprintf(fp, "graph G {\n");
    for (i=0; i<m_c_idx.size(); i++) {
      p = m_c_idx[i];


      for (j=0; j<m_E[p].size(); j++) {
        if (j==0) {
          fprintf(fp, "  c%i [shape=box];\n", m_E[p][j].m_id[0]);
        }

        fprintf(fp, "  c%i -- v%i", m_E[p][j].m_id[0], m_E[p][j].m_id[1]);
        if (m_E[p][j].m_J > 0) {
          fprintf(fp, " [style=dotted,label=\"%0.2f\"]", (float)m_E[p][j].m_d[m_t][0]);
        }
        else {
          fprintf(fp, " [label=\"%0.2f\"]", (float)m_E[p][j].m_d[m_t][0]);
        }
        fprintf(fp, ";\n");
      }
    }
    fprintf(fp, "}\n");
  }

  int print_dot_w_fn(char *fn) {
    FILE *fp;
    fp=fopen(fn, "w");
    if (!fp) { return -1; }
    print_dot_w(fp);
    fclose(fp);
    return 0;
  }

  int tick(void) {
    int32_t i, j, n;
    int32_t t_cur, t_nxt, J;
    int32_t v_idx, c_idx, p;
    int32_t _q, _s;
    double prod_plus, prod_minus, prod_all;
    double prod_u, prod_s, prod_0;
    double d;

    t_cur = m_t;
    t_nxt = 1-m_t;

    // first update variable 'probabilities'
    //
    for (v_idx=0; v_idx < m_v_idx.size(); v_idx++) {
      p = m_v_idx[v_idx];

      prod_plus = 1.0;
      prod_minus = 1.0;
      prod_all = 1.0;
      for (i=0; i<m_E[p].size(); i++) {

        _q = m_E[p][i].m_idx[1];
        _s = m_E[p][i].m_bp_offset;

        d = m_E[_q][_s].m_d[t_cur][0];
        prod_all *= (1.0 - d);
        if (m_E[p][i].m_J < 0) {
          prod_plus *= (1.0 - d);
        }
        else {
          prod_minus *= (1.0 - d);
        }
      }

      for (i=0; i<m_E[p].size(); i++) {

        _q = m_E[p][i].m_idx[1];
        _s = m_E[p][i].m_bp_offset;

        d = (1.0 - m_E[_q][_s].m_d[t_cur][0]);
        if (d<m_eps) { d=1.0; }

        //>>UNCHECKED!!!!!!!
        if ( m_E[p][i].m_J == 1 ) {
          m_E[p][i].m_d[t_nxt][0] = (1.0 - prod_plus) * (prod_minus / d);
          m_E[p][i].m_d[t_nxt][1] = (1.0 - (prod_minus / d) ) * prod_plus;
        }
        else {
          m_E[p][i].m_d[t_nxt][0] = (1.0 - prod_minus) * (prod_plus / d);
          m_E[p][i].m_d[t_nxt][1] = (1.0 - (prod_plus / d) ) * prod_minus;
        }
        m_E[p][i].m_d[t_nxt][2] = prod_all / d;
        //<<UNCHECKED!!!!!!!

      }

    }

    // next, using the above, update the clause 'probabilities'
    //
    for (c_idx=0; c_idx < m_c_idx.size(); c_idx++) {
      p = m_c_idx[c_idx];

      prod_u = 1.0;
      prod_all = 1.0;
      for (i=0; i<m_E[p].size(); i++) {

        _q = m_E[p][i].m_idx[1];
        _s = m_E[p][i].m_bp_offset;

        prod_u *= m_E[_q][_s].m_d[t_nxt][0];
        prod_all *= ( m_E[_q][_s].m_d[t_nxt][0] +
                      m_E[_q][_s].m_d[t_nxt][1] +
                      m_E[_q][_s].m_d[t_nxt][2] );
        
      }

      if (prod_all < m_eps) { prod_all = 1.0; }

      for (i=0; i<m_E[p].size(); i++) {
        _q = m_E[p][i].m_idx[1];
        _s = m_E[p][i].m_bp_offset;

        d = m_E[_q][_s].m_d[t_nxt][0];
        if (d < m_eps) {
          m_E[p][i].m_d[t_nxt][0] = 1.0;
          continue;
        }

        d = ( m_E[_q][_s].m_d[t_nxt][0] +
              m_E[_q][_s].m_d[t_nxt][1] +
              m_E[_q][_s].m_d[t_nxt][2] ) / d;
        m_E[p][i].m_d[t_nxt][0] = (prod_u / prod_all) * d;
      }
    }

    m_t = t_nxt;
  }

  void print_sp(void) {
  }

  void print(void) {
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
  int32_t ch, i, j;
  std::string buf;
  std::vector<int32_t> int_line;
  int32_t _start, _n;

  int32_t _v_n, _v_s, _v_idx;
  int32_t _c_n, _c_s, _c_idx;

  int32_t varname;

  int32_t clause_base_idx, var_base_idx;
  int32_t c_idx, c_n, v_idx, v_n;
  int32_t clause_idx, clause_len;
  int32_t var_idx;
  int32_t _idx;

  sp_edge_t _edge;

  std::vector<int32_t> v_nvar;

  while (!feof(fp)) {
    ch = fgetc(fp);
    if ((ch==EOF) || (ch=='\n')) {

      if ( (buf.size()>0) &&
           (((buf[0] == '-') ||
            ((buf[0] >= '0') && (buf[0] <= '9')))) ) {
        int_line.clear();
        parse_int_line(buf, int_line);

        _start = sp.m_raw_clause.size();
        _n = 0;
        sp.m_raw_clause_info.push_back(_start);
        for (i=0; i<int_line.size(); i++) {
          if (int_line[i] == 0) { break; }
          sp.m_raw_clause.push_back(int_line[i]);

          if (abs(int_line[i]) >= sp.m_nvar) {
            sp.m_nvar = abs(int_line[i]);
          }
          _n++;
        }
        sp.m_raw_clause_info.push_back(_n);
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

    _start = sp.m_raw_clause.size();
    _n = 0;
    sp.m_raw_clause_info.push_back(_start);
    for (i=0; i<int_line.size(); i++) {
      if (int_line[i] == 0) { break; }
      sp.m_raw_clause.push_back(int_line[i]);

      if (abs(int_line[i]) >= sp.m_nvar) {
        sp.m_nvar = abs(int_line[i]);
      }
    }
    sp.m_raw_clause_info.push_back(_n);
    sp.m_nclause++;
  }

  //---
  //---

  // parsing is done, do edge creation
  //

  sp.m_c_idx.resize( sp.m_nclause );
  sp.m_v_idx.resize( sp.m_nvar );
  sp.m_E.resize( sp.m_nclause + sp.m_nvar );

  clause_base_idx = 0;
  var_base_idx = sp.m_nclause;

  _edge.m_d[0][0] = 0.0;
  _edge.m_d[0][1] = 0.0;
  _edge.m_d[0][2] = 0.0;

  _edge.m_d[1][0] = 0.0;
  _edge.m_d[1][1] = 0.0;
  _edge.m_d[1][2] = 0.0;

  for (i=0; i<sp.m_c_idx.size(); i++) {
    sp.m_c_idx[i] = clause_base_idx + i;
  }
  for (i=0; i<sp.m_v_idx.size(); i++) {
    sp.m_v_idx[i] = var_base_idx + i;
  }

  for (i=0; i<sp.m_raw_clause_info.size(); i+=2) {

    clause_idx = i/2;
    c_idx = sp.m_raw_clause_info[i];
    c_n   = sp.m_raw_clause_info[i+1];

    sp.m_c_idx[clause_idx] = clause_base_idx + clause_idx;

    for (j=0; j<c_n; j++) {
      varname = abs(sp.m_raw_clause[c_idx + j]);
      var_idx = varname-1;

      _edge.m_type = (int)'c';
      _edge.m_id[0]= clause_idx;
      _edge.m_id[1]= varname;
      _edge.m_idx[0] = clause_base_idx + clause_idx;
      _edge.m_idx[1] = var_base_idx + var_idx;
      _edge.m_J = ( (sp.m_raw_clause[c_idx + j] < 0) ? 1 : -1 );
      _edge.m_d[0][0] = _drand();
      sp.m_E[clause_base_idx + clause_idx].push_back(_edge);

      _edge.m_type = (int)'v';
      _edge.m_id[0] = varname;
      _edge.m_id[1] = clause_idx;
      _edge.m_idx[0] = var_base_idx + var_idx;
      _edge.m_idx[1] = clause_base_idx + clause_idx;
      _edge.m_J = ( (sp.m_raw_clause[c_idx + j] < 0) ? 1 : -1 );
      _edge.m_d[0][0] = _drand();
      _edge.m_d[0][1] = _drand();
      _edge.m_d[0][2] = _drand();
      _edge.m_d[1][0] = 0.0;
      _edge.m_d[1][1] = 0.0;
      _edge.m_d[1][2] = 0.0;
      sp.m_E[var_base_idx + var_idx].push_back(_edge);

      _c_idx = sp.m_E[clause_base_idx + clause_idx].size()-1;
      _v_idx = sp.m_E[var_base_idx + var_idx].size()-1;
      sp.m_E[clause_base_idx + clause_idx][_c_idx].m_bp_offset = _v_idx;
      sp.m_E[var_base_idx + var_idx][_v_idx].m_bp_offset = _c_idx;
    }

  }

  return 0;
}

int main(int argc, char **argv) {
  int i, ch, opt_idx, r;
  FILE *ifp = stdin;
  std::string ifn;
  std::string ofn;

  char _ofn[1024];

  std::vector<int32_t> nei;

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

  r = sp.consistency();
  //printf("# got %i\n", r);

  //sp.print_E();
  //sp.print_dimacs(stdout);

  for (i=0; i<10; i++) {
    snprintf(_ofn, 32, "./tmp/t%02i.dot", i);
    sp.print_dot_w_fn(_ofn);
    sp.tick();
  }
  //sp.print_dot_w(stdout);
  //sp.print_dot(stdout);

  //sp.tick();


  //sp.print_dimacs(stdout);
  //sp.print_E();

  //sp.init_sp();

  //sp.consistency();

  //sp.print();
  //sp.print_sp();



}
