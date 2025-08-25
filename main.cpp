#include "MyLib/MyTimer.h"
#include "MyLib/ParaReader.h"
#include "graph/Graph.h"
#include <iostream>
#include <set>
#include <stdio.h>
#include <random>

struct Param {
    char graph_file[200];
    char update_file[200];
    int w_type;
    float a;
    float b;
};

struct Param parseArgs(int nargs, char **args) {
    Param rtn;
    int cnt = 1;
    bool failed = false;
    char *arg;
    int i;
    char para[10];
    char graph_file[200] = "";
    char update_file[200] = "";
    int w_type = 1;
    float a = 0.5;
    float b = 0.5;
//    printf("The input parameters are:\n\n");
    while (cnt < nargs && !failed) {
        arg = args[cnt++];
        if (cnt == nargs) {
            failed = true;
            break;
        }
        i = getNextChar(arg);
        if (arg[i] != '-') {
            failed = true;
            break;
        }
        getNextWord(arg + i + 1, para);
        arg = args[cnt++];
        if (strcmp(para, "graph") == 0) {
            getNextWord(arg, graph_file);
//            printf("%s\n", graph_file);
        } else if (strcmp(para, "update") == 0) {
            getNextWord(arg, update_file);
//            printf("%s\n", update_file);
        } else if (strcmp(para, "wtype") == 0) {
            w_type = atof(arg);
            printf("weight type : %d\n", w_type);
        } else if (strcmp(para, "a") == 0) {
            a = atof(arg);
            printf("a : %lf\n", a);
        } else if (strcmp(para, "b") == 0) {
            b = atof(arg);
            printf("b : %lf\n", b);
        } else {
            failed = true;
            printf("Unknown option -%s!\n\n", para);
        }
    }

    /*****************************************************************************/
    strcpy(rtn.graph_file, graph_file);
    strcpy(rtn.update_file, update_file);
    rtn.a = a;
    rtn.b = b;
    rtn.w_type = w_type;
    return rtn;
}

void usage() {
    printf("Usage:\n");
    printf("dagp -graph [graph file] -update [update file] ");
}

int main(int argc, char **argv) {
    if (argc <= 2) {
        usage();
        return 0;
    }
    printf("Start to parse the arguments\n");

    Param para = parseArgs(argc, argv);

    printf("Arguments parsed.\n");

    printf("---------------------------------------------------------------------\n");

    FILE *f = fopen(para.graph_file, "rb");

    if (f == NULL) {
        printf("graph file not found.\n");
        exit(1);
    }

    unsigned int n, m;
    fread(&n, 1, sizeof(int), f);
    fread(&m, 1, sizeof(int), f);
    int *edges = (int *) malloc(sizeof(int) * m);
    fread(edges, m, sizeof(int), f);
    fclose(f);

    MyVector<dynagp::Vertex *> _vList;
    _vList.reserve(n);

    for (int i = 0; i < n; i++) {
        dynagp::Vertex *newVertex = new dynagp::Vertex(i + 1);
        _vList.push_back(newVertex);
    }
    Graph graph(_vList);

    _vList.release_space();

    for (int i = 0, j = 0; i < m; i += 2, ++j) {
        graph.insertEdge(edges[i], edges[i + 1]);
    }
    printf("Graph input with %d vertices and %d edges.\n", n,
           m / 2);
    printf("---------------------------------------------------------------------\n");

    f = fopen(para.update_file, "r");
    vector<pair<int, pair<int, int>>> updates;
    updates.reserve(9 * (m / 2));
    int c = 0, d = 0, e = 0;
    while (fscanf(f, "%d%d%d", &c, &d, &e) != EOF) {
        updates.emplace_back(make_pair(c, make_pair(d, e)));
    }
    fclose(f);

    printf("%ld Graph update file loaded.\n", updates.size());
    printf("---------------------------------------------------------------------\n");
    printf("initialization test:\n");
    //initialization test
    std::vector<double> x(n);
    std::vector<double> r(n);
    double start = 0;
    double end = 0;
    double total_time = 0;
    double delta = 0;
    double eps = 0;
    int L = (int) (log2(n));
    // test delta from 1/100 to 1/1000
     for (int i = 100; i <= 1000; i = i + 100) {
         delta = 1.0 / i;
         eps = 0.1 * 0.1 * delta / (L + 1) / 2;
         std::random_device rd;
         std::mt19937 gen(rd());
         std::uniform_int_distribution<int> dis(1, (int) (L / delta));
         int group_num = dis(gen);
         std::vector<int> group_size(group_num);
         std::vector<double> group_value(group_num);
         int left = n;
         for (int j = 0; j < group_num - 1; j++) {
             std::uniform_int_distribution<int> dis2(1, left - (group_num - j - 1));
             group_size[j] = dis2(gen);
             left -= group_size[j];
         }
         group_size[group_num - 1] = left;

         double left_v = 1.0;
         for (int j = 0; j < group_num; j++) {
             std::uniform_real_distribution<double> realDist2(0.0, left_v);
             group_value[j] = realDist2(gen);
             left_v = left_v - group_value[j];
             group_value[j] = group_value[j] / group_size[j];
         }

         start = getCurrentTime();
         int z = 0;
         double value = 0;
         for (int j = 0; j < group_num; j++) {
             for (int k = 0; k < group_size[j]; k++) {
                 x[z] = group_value[j];
                 z = z + 1;
                 value = value + group_value[j];
             }
         }
         end = getCurrentTime();
         total_time = end - start;
         printf("Initialization time %.5f: %.9lf\t",
                delta, total_time);
         z = 0;
         std::uniform_real_distribution<double> realDist(0.0, 1.0);
         start = getCurrentTime();
         for (int j = 0; j < group_num; ++j) {
             double gv = group_value[j];
             int gs = group_size[j];

             if (gv <= 0.0) continue;

             double p = gv / eps;
             if (p >= 1.0) {
                 std::fill(r.begin() + z, r.begin() + z + gs, gv);
                 z += gs;
             } else {
                 double log2_p_inv = 1.0 / std::log2(p);  // precompute
                 int o = 0;
                 while (o < gs) {
                     double u = realDist(gen);
                     int y = static_cast<int>(std::ceil(std::log2(u) * log2_p_inv));
                     o += y;
                     if (o >= gs) break;
                     r[z + o] = eps;
                 }
                 z += gs;
             }
         }

         end = getCurrentTime();
         total_time = end - start;
         printf("%.9lf\n",
                total_time);
     }
     printf("---------------------------------------------------------------------\n");
//query test
    printf("Query test:\n");
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> realDist(0.0, 1.0);
    std::uniform_int_distribution<int> dist(0, n-1);
    float a = para.a;
    float b = para.b;
    // float a = 0;
    // float b = 1;

    auto *x_ = new double[n];
    for (int i = 0; i < n; i++) {
        x_[i] = 1.0/n;
        // x_[i] = realDist(gen);
    }

    // normalize
    // double sum = 0;
    // for (int i = 0; i < n; i++) {
    //     sum = sum + x_[i];
    // }
    // for (int i = 0; i < n; i++) {
    //     x_[i] = 1.0 / sum;
    // }

    //one hot
    // for (int i = 0; i < n; i++) {
    //     x_[i] = 0;
    // }
    // x_[dist(gen)] = 1.0;

    L = (int) (log2(n));
    delta = 1.0/n;
    // L = L * 2;
    auto *w = new double[L + 1];
    if (para.w_type == 0) {
        for (int i = 0; i < L; i++) {
            w[i] = 0;
        }
        w[L] = 1;
    } else if (para.w_type == 1) {
        float alpha = 0.2;
        w[0] = alpha;
        for (int i = 1; i < L + 1; i++) {
            w[i] = w[i - 1] * (1 - alpha);
        }
    } else if (para.w_type == 2) {
        float lambda = 1;
        float beta = 0.85 / lambda;
        w[0] = 1 - beta;
        for (int i = 1; i < L + 1; i++) {
            w[i] = w[i - 1] * beta;
        }
    } else if (para.w_type == 3) {
        float t = 4;
        w[0] = exp(-t);
        for (int i = 1; i < L + 1; i++) {
            w[i] = w[i - 1] * t / i;
        }
    }
    // printf("L:%d\n", L);
    eps = 0.1 * 0.1 * delta / (L + 1) / 2;
    // eps = 0.1 * 0.1 * delta / L / L;
    // printf("eps:%.9f\n", eps);
    graph.query(x_, a, b, L, w, eps, 'S');


//update test
    a = 0.5;
    b = 0.5;
    double total_update = 0;
    graph.construct_DS(a, b);
    graph.DAGP_initialize();
    for (auto &update: updates) {
        start = getCurrentTime();
        if (update.first == 1) {
//            graph.naive_insert(update.second.first, update.second.second, a, b);
            graph.DAGP_insert(update.second.first, update.second.second, a, b);
        } else {
//            graph.naive_delete(update.second.first, update.second.second, a, b);
            graph.DAGP_delete(update.second.first, update.second.second, a, b);
        }
        end = getCurrentTime();
        total_update += end - start;
    }

    printf("---------------------------------------------------------------------\n");
    printf("total update time: %.9lf\t",
           total_update);

}


