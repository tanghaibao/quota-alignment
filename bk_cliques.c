/* Author: Ram Samudrala (me@ram.org) April 19, 1995
 *         extend() was implemented from the Bron and Kerbosch
 *         algorithm (Communications of the ACM 16:575-577, 1973) by
 *         Jean-Francios Gibrat.
 *
 * Modified by Haibao Tang (bao@uga.edu) April 26, 2007
 * 1. use dynamic memory alloc
 * 2. use edge list since most graphs we deal with are sparse
 *
 * Input file format is the following:
 * (assume number of nodes N, size-K cliques to report)
 *
 * N K
 * v1 v2
 * v3 v4
 * ...
 * pairs of nodes numbered between 0 and N-1
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AllocArray(pt, size) (pt = malloc(sizeof(*pt) * (size)))

static int num_nodes, K, nel, *all, *compsub;
static char **connected;

static void init_memory();
static void read_graph();
static int extend(int old[], int ne, int ce);


int main(int argc, char *argv[])
{
    read_graph();

    for (nel = 1; nel <= num_nodes; nel++)
        all[nel] = nel;
    nel = 0;

    extend(all, 0, num_nodes);

    return 0;

}


void init_memory()
/*  initialize memory  */
{
    int i;
    AllocArray(all, num_nodes);
    AllocArray(compsub, num_nodes);
    AllocArray(connected, num_nodes+1);
    for (i = 0; i <= num_nodes; i++)
        AllocArray(connected[i], num_nodes+1);

}


void read_graph()
{
    int i, j, a, b;

    /* first row contains N (number of nodes) and K (minimum size clique to report) */
    scanf("%d%d", &num_nodes, &K);
    fprintf(stderr, "Number of nodes in the graph: %d\n", num_nodes);

    init_memory();

    for (i = 0; i <= num_nodes; i++)
    {
        for (j = 0; j <= num_nodes; j++)
            connected[i][j] = '0';
        connected[i][i] = '1';
    }

    while (scanf("%d%d", &a, &b)==2)
    {

        a++;
        b++;

        connected[a][b] = '1';
        connected[b][a] = '1';
    }

#if 0
    for (i = 1; i <= num_nodes; i++)
    {
        for (j = 1; j <= num_nodes; j++)
            fprintf(stdout, "%c", connected[i][j]);
        fprintf(stdout, "\n");
    }
    exit(1);
#endif

}



int extend(int *old, int ne, int ce)
/*
 * Author: Jean-Francois Gibrat 13/07/94
 * Modified: Haibao Tang <bao@uga.edu> April. 27, 2007
 *
 * This subroutine finds the maximal complete subgraphs (cliques) of a
 * given undirected graph. The algorithm consists in a recursively defined
 * extension operator (extend) applied to the three following sets:
 * set compsub: set to be extended by a new point or shrunk by one point as
 *              the algorithm travels along a branch of the backtracking tree.
 * set cand:    set of all vertices connected to the current vertex in compsub
 *              i.e., vertices which will be used as an extension of the
 *              present configuration of compsub.
 * set not:     set of all points which have already served as an extension and
 *              are now excluded.
 *
 * The basic mechanism consists in the following steps:
 * Step 1:  selection of a candidate.
 * Step 2:  creating new sets cand and not from the old sets by removing
 *       :  all points not connected to the selected candidate,
 *       :  keeping the old sets intact.
 * Step 3:  adding the selected candidate to compsub.
 * Step 4:  calling the extension operator to operate on the sets just formed.
 * Step 5:  upon return the selected candidate is transfered from compsub
 *       :  to the old set not.
 *
 * The set compsub is a stack external to the subroutine.
 * The sets not and cand are stored in the local array old,
 * not from positions 1 to ne and cand from positions ne+1 to ce.
 * nvert x nvert array [connected] indicates which vertices are connected
 * to which vertices. This array is external to the subroutine.
 *
 * Reference: Bron C. and Kerbosch J., Commun.
 * A.C.M., vol. 16, pp 575-577, 1973.
 *
 * Comments by Haibao Tang:
 * The variable names were taken directly from the published pseudocode,
 * indeed this may cause confusion. Nemonics are the following:
 * ne: not end; ce: cand(didate) end; minnod: min # of disconnections;
 * sel: selected;
 *
 * the not and cand sets are stored in a 1D array with this layout
 * | not |  cand  |
 * 0 ... ne ... ce
 *
 *
 */

{
    int s = 0;
    int nod = 0;
    int pos = 0;
    int fixp = 0;
    int minnod = ce;
    int knod, newne, newce, count, p, sel;
    size_t i, j;
    int *_new;
    AllocArray(_new, num_nodes);

    /*
     * Look for the vertex that is connected to the maximum number of vertices
     * i.e., the minimum number of disconnections (Step 1)
     *
     */
    for (i = 1; (i <= ce) && (minnod != 0); i++)
    {
        p = old[i];
        count = 0;
        /* Count disconnections */
        for (j = ne + 1; (j <= ce) && (count < minnod); j++)
        {
            if (connected[p][old[j]] == '0')
            {
                count++;
                pos = j;  /* save position of potential candidate */
            }
        }
        /* Test whether vertex at position i in the stack has a minimum
         * number of diconnections
         */
        if (count < minnod)
        {
            fixp = p;
            minnod = count;
            if (i <= ne)
                s = pos;
            else
            {
                s = i;
                nod = 1;  /* Number of disconnections preincreased by one */
            }
        }
    }

    /*
     * Backtrack Cycle: the loop is executed as many times as
     * there are vertices not connected to the selected point
     */

    for (knod = minnod + nod; knod >= 1; knod--)
    {
        /*
         * Interchange
         */
        p = old[s];
        old[s] = old[ne+1];
        old[ne+1] = p;
        sel = old[ne+1];
        /*
         * Fill new set not (Step 2)
         */
        newne = 0;
        for (i = 1; i <= ne; i++)
            if (connected[sel][old[i]] == '1')
                _new[++newne] = old[i];
        /*
         * Fill new set cand (Step 2)
         */
        newce = newne;
        for (i = ne + 2; i <= ce; i++)
            if (connected[sel][old[i]] == '1')
                _new[++newce] = old[i];

        /*
         * The selected vertex is added to compsub (Step 3)
         */
        compsub[++nel] = sel; /* Note index nel is external to the subroutine */

        if (newce == 0)
        {
            /* The new clique is stored in compsub. */
            if (nel > K)
            {
                for (i = 1; i <= nel; i++) printf("%d ", compsub[i]-1);
                printf("\n");
            }
        }
        /*
         * Otherwise calls the extension operator to operate on the sets
         * just formed  (Step 4)
         */
        else
        {
            if (newne < newce)
                extend(_new, newne, newce);
        }

        /*
         * Upon return removal of the selected candidate from compsub
         * and its addition to the old set not (Step 5)
         *
         * Remove from compsub */
        nel--;
        /* Add to not */
        ne++;

        /*
         * Select a candidate disconnected to the fixed point (Step 1)
         */
        if ( knod > 1)
        {
            s = ne + 1;
            while (connected[fixp][old[s]] == '1')
                s++;
        }
    } /* End of the loop with index knod */

    free(_new);

    return 1;

}

