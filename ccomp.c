
//
// CCOMP III
// receptor/complex comparison program
// by piotr rotkiewicz, piotr /at/ pirx /dot/ com, august 2006
//
// Version 3.70
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FLAG_BACKBONE 1
#define FLAG_SIDECHAIN 2

typedef struct _atom_struct {
    int anum;
    int rnum;
    int flags;
    char aa;
    char rid[7];
    char aname[6];
    char rname[5];
    float x, y, z;
    char het;
    char chain;
    float dat1, dat2;
    struct _atom_struct* next;
} atom_struct;

typedef struct _mol_struct {
    atom_struct* atoms;
    int natom;
    float res;
    struct _mol_struct* next;
} mol_struct;

short int** cmap;

char AA_NAMES[21][4] = { "GLY", "ALA", "SER", "CYS", "VAL",
    "THR", "ILE", "PRO", "MET", "ASP",
    "ASN", "LEU", "LYS", "GLU", "GLN",
    "ARG", "HIS", "PHE", "TYR", "TRP",
    "UNK" };

char SHORT_AA_NAMES[22] = { "GASCVTIPMDNLKEQRHFYWX" };

int nheavy[21] = { 0, 1, 2, 2, 3, 3, 4, 3, 4, 4, 4, 4, 5, 5, 5, 7, 6, 7, 8, 10, 0 };

char* heavy_atoms[200] = {
    /* GLY */ NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    /* ALA */ "CB ", NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    /* SER */ "CB ", "OG ", NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    /* CYS */ "CB ", "SG ", NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    /* VAL */ "CB ", "CG1", "CG2", NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    /* THR */ "CB ", "OG1", "CG2", NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    /* ILE */ "CB ", "CG1", "CG2", "CD1", NULL, NULL, NULL, NULL, NULL, NULL,
    /* PRO */ "CB ", "CG ", "CD ", NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    /* MET */ "CB ", "CG ", "SD ", "CE ", NULL, NULL, NULL, NULL, NULL, NULL,
    /* ASP */ "CB ", "CG ", "OD1", "OD2", NULL, NULL, NULL, NULL, NULL, NULL,
    /* ASN */ "CB ", "CG ", "OD1", "ND2", NULL, NULL, NULL, NULL, NULL, NULL,
    /* LEU */ "CB ", "CG ", "CD1", "CD2", NULL, NULL, NULL, NULL, NULL, NULL,
    /* LYS */ "CB ", "CG ", "CD ", "CE ", "NZ ", NULL, NULL, NULL, NULL, NULL,
    /* GLU */ "CB ", "CG ", "CD ", "OE1", "OE2", NULL, NULL, NULL, NULL, NULL,
    /* GLN */ "CB ", "CG ", "CD ", "OE1", "NE2", NULL, NULL, NULL, NULL, NULL,
    /* ARG */ "CB ", "CG ", "CD ", "NE ", "CZ ", "NH1", "NH2", NULL, NULL, NULL,
    /* HIS */ "CB ", "CG ", "ND1", "CD2", "CE1", "NE2", NULL, NULL, NULL, NULL,
    /* PHE */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ", NULL, NULL, NULL,
    /* TYR */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ", "OH ", NULL, NULL,
    /* TRP */ "CB ", "CG ", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"
};

int AA_NUMS[256];

char setseq(char* aaname)
{
    int i;

    for (i = 0; i < 21; i++)
        if ((aaname[0] == AA_NAMES[i][0]) && (aaname[1] == AA_NAMES[i][1]) && (aaname[2] == AA_NAMES[i][2]))
            break;
    if (i == 21) {
        if (!strcmp(aaname, "GLX"))
            return 'E';
        if (!strcmp(aaname, "ASX"))
            return 'D';
        if (!strcmp(aaname, "MSE"))
            return 'M';
        if (!strcmp(aaname, "SEP"))
            return 'S';
        if (!strcmp(aaname, "TPO"))
            return 'T';
        if (!strcmp(aaname, "PTR"))
            return 'Y';
        i--;
    }

    return SHORT_AA_NAMES[i];
}

mol_struct* read_pdb(char* name)
{
    FILE* inp;
    int natom, nmol;
    float res;
    char *ptr, buf[1000];
    atom_struct *atom, *lastatom;
    mol_struct *mol, *lastmol, *firstmol;

    //  printf("Reading %s...\n", name);

    natom = nmol = 0;
    atom = lastatom = NULL;
    mol = lastmol = firstmol = NULL;

    inp = fopen(name, "r");
    if (!inp)
        return NULL;

    res = -1.0;
    while (!feof(inp)) {
        if (fgets(buf, 1000, inp) == buf) {
            if ((!strncmp(buf, "REMARK", 6) && strstr(buf, "RESOLUTION."))) {
                ptr = strstr(buf, "RESOLUTION.");
                if (!sscanf(ptr + 11, "%f", &res))
                    res = -1.0;
            }
            if ((!strncmp(buf, "ATOM", 4) || !strncmp(buf, "HETATM", 6)) && !(buf[17] == 'H' && buf[18] == 'O' && buf[19] == 'H')) {
                //printf("%c%c\n", buf[12], buf[13]);
                atom = (atom_struct*)calloc(sizeof(atom_struct), 1);
                if (buf[13] == 'N' && buf[14] == ' ')
                    atom->flags |= FLAG_BACKBONE;
                else if (buf[13] == 'C' && buf[14] == ' ')
                    atom->flags |= FLAG_BACKBONE;
                else if (buf[13] == 'O' && buf[14] == ' ')
                    atom->flags |= FLAG_BACKBONE;
                else
                    atom->flags |= FLAG_SIDECHAIN;
                sscanf(&buf[6], "%d", &atom->anum);
                sscanf(&buf[22], "%d", &atom->rnum);
                strncpy(atom->rid, &buf[22], 5);
                atom->rid[6] = 0;
                strncpy(atom->aname, &buf[12], 4);
                atom->aname[5] = 0;
                strncpy(atom->rname, &buf[17], 3);
                atom->rname[4] = 0;
                sscanf(&buf[30], "%f", &atom->x);
                sscanf(&buf[38], "%f", &atom->y);
                sscanf(&buf[46], "%f", &atom->z);
                sscanf(&buf[54], "%f", &atom->dat1);
                sscanf(&buf[60], "%f", &atom->dat2);
                atom->chain = buf[21];
                atom->aa = setseq(atom->rname);
                if (!strncmp(buf, "HETATM", 6))
                    atom->het = 1;
                natom++;
                if (!mol) {
                    mol = (mol_struct*)calloc(sizeof(mol_struct), 1);
                    mol->atoms = lastatom = atom;
                }
                else {
                    lastatom->next = atom;
                    lastatom = atom;
                }
            }
            if ((!strncmp(buf, "TER", 3) || !strncmp(buf, "END", 3)) && natom) {
                ++nmol;
                //   printf("Molecule %2d includes %4d atoms\n", ++nmol, natom);
                mol->natom = natom;
                mol->res = res;
                natom = 0;
                if (lastmol) {
                    lastmol->next = mol;
                }
                else
                    firstmol = mol;
                lastmol = mol;
                mol = NULL;
            }
        }
    }
    if (natom) {
        ++nmol;
        //  printf("Molecule %2d includes %4d atoms\n", ++nmol, natom);
        mol->natom = natom;
        mol->res = res;
        if (lastmol) {
            lastmol->next = mol;
        }
        else
            firstmol = mol;
    }
    fclose(inp);

    return firstmol;
}

float superimpose(float** coords1, float** coords2, int npoints, float** rot, float* trans, float* cm)
{
    float mat_s[3][3], mat_a[3][3], mat_b[3][3], mat_g[3][3];
    float mat_u[3][3], tmp_mat[3][3];
    float val, d, alpha, beta, gamma, x, y, z;
    float cx1, cy1, cz1, cx2, cy2, cz2, tmpx, tmpy, tmpz;
    int i, j, k, n;

    cx1 = cy1 = cz1 = cx2 = cy2 = cz2 = 0.;

    for (i = 0; i < npoints; i++) {
        cx1 += coords1[i][0];
        cy1 += coords1[i][1];
        cz1 += coords1[i][2];
        cx2 += coords2[i][0];
        cy2 += coords2[i][1];
        cz2 += coords2[i][2];
    }

    cx1 /= (float)npoints;
    cy1 /= (float)npoints;
    cz1 /= (float)npoints;

    cx2 /= (float)npoints;
    cy2 /= (float)npoints;
    cz2 /= (float)npoints;

    for (i = 0; i < npoints; i++) {
        coords1[i][0] -= cx1;
        coords1[i][1] -= cy1;
        coords1[i][2] -= cz1;
        coords2[i][0] -= cx2;
        coords2[i][1] -= cy2;
        coords2[i][2] -= cz2;
    }

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            if (i == j)
                mat_s[i][j] = mat_a[i][j] = mat_b[i][j] = mat_g[i][j] = 1.0;
            else
                mat_s[i][j] = mat_a[i][j] = mat_b[i][j] = mat_g[i][j] = 0.0;
            mat_u[i][j] = 0.;
        }

    for (n = 0; n < npoints; n++) {
        mat_u[0][0] += coords1[n][0] * coords2[n][0];
        mat_u[0][1] += coords1[n][0] * coords2[n][1];
        mat_u[0][2] += coords1[n][0] * coords2[n][2];
        mat_u[1][0] += coords1[n][1] * coords2[n][0];
        mat_u[1][1] += coords1[n][1] * coords2[n][1];
        mat_u[1][2] += coords1[n][1] * coords2[n][2];
        mat_u[2][0] += coords1[n][2] * coords2[n][0];
        mat_u[2][1] += coords1[n][2] * coords2[n][1];
        mat_u[2][2] += coords1[n][2] * coords2[n][2];
    }

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            tmp_mat[i][j] = 0.;

    do {
        d = mat_u[2][1] - mat_u[1][2];
        if (d == 0)
            alpha = 0;
        else
            alpha = atan(d / (mat_u[1][1] + mat_u[2][2]));
        if (cos(alpha) * (mat_u[1][1] + mat_u[2][2]) + sin(alpha) * (mat_u[2][1] - mat_u[1][2]) < 0.0)
            alpha += M_PI;
        mat_a[1][1] = mat_a[2][2] = cos(alpha);
        mat_a[2][1] = sin(alpha);
        mat_a[1][2] = -mat_a[2][1];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    tmp_mat[i][j] += mat_u[i][k] * mat_a[j][k];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                mat_u[i][j] = tmp_mat[i][j];
                tmp_mat[i][j] = 0.;
            }
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    tmp_mat[i][j] += mat_a[i][k] * mat_s[k][j];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                mat_s[i][j] = tmp_mat[i][j];
                tmp_mat[i][j] = 0.;
            }
        d = mat_u[0][2] - mat_u[2][0];
        if (d == 0)
            beta = 0;
        else
            beta = atan(d / (mat_u[0][0] + mat_u[2][2]));
        if (cos(beta) * (mat_u[0][0] + mat_u[2][2]) + sin(beta) * (mat_u[0][2] - mat_u[2][0]) < 0.0)
            beta += M_PI;
        mat_b[0][0] = mat_b[2][2] = cos(beta);
        mat_b[0][2] = sin(beta);
        mat_b[2][0] = -mat_b[0][2];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    tmp_mat[i][j] += mat_u[i][k] * mat_b[j][k];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                mat_u[i][j] = tmp_mat[i][j];
                tmp_mat[i][j] = 0.;
            }
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    tmp_mat[i][j] += mat_b[i][k] * mat_s[k][j];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                mat_s[i][j] = tmp_mat[i][j];
                tmp_mat[i][j] = 0.;
            }
        d = mat_u[1][0] - mat_u[0][1];
        if (d == 0)
            gamma = 0;
        else
            gamma = atan(d / (mat_u[0][0] + mat_u[1][1]));
        if (cos(gamma) * (mat_u[0][0] + mat_u[1][1]) + sin(gamma) * (mat_u[1][0] - mat_u[0][1]) < 0.0)
            gamma += M_PI;
        mat_g[0][0] = mat_g[1][1] = cos(gamma);
        mat_g[1][0] = sin(gamma);
        mat_g[0][1] = -mat_g[1][0];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    tmp_mat[i][j] += mat_u[i][k] * mat_g[j][k];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                mat_u[i][j] = tmp_mat[i][j];
                tmp_mat[i][j] = 0.;
            }
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    tmp_mat[i][j] += mat_g[i][k] * mat_s[k][j];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                mat_s[i][j] = tmp_mat[i][j];
                tmp_mat[i][j] = 0.;
            }
        val = fabs(alpha) + fabs(beta) + fabs(gamma);
    } while (val > 1e-3);

    val = 0.;
    for (i = 0; i < npoints; i++) {
        x = coords2[i][0];
        y = coords2[i][1];
        z = coords2[i][2];
        tmpx = x * mat_s[0][0] + y * mat_s[0][1] + z * mat_s[0][2];
        tmpy = x * mat_s[1][0] + y * mat_s[1][1] + z * mat_s[1][2];
        tmpz = x * mat_s[2][0] + y * mat_s[2][1] + z * mat_s[2][2];
        x = coords1[i][0] - tmpx;
        y = coords1[i][1] - tmpy;
        z = coords1[i][2] - tmpz;
        val += x * x + y * y + z * z;
    }

    trans[0] = cx1;
    trans[1] = cy1;
    trans[2] = cz1;

    cm[0] = cx2;
    cm[1] = cy2;
    cm[2] = cz2;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            rot[i][j] = mat_s[i][j];

    for (i = 0; i < npoints; i++) {
        coords1[i][0] += cx1;
        coords1[i][1] += cy1;
        coords1[i][2] += cz1;
        coords2[i][0] += cx2;
        coords2[i][1] += cy2;
        coords2[i][2] += cz2;
    }

    return sqrt(val / (float)npoints);
}

float **SCORE, **VAL;
char** MOV;

float GAP_OPEN_PENALTY = -1.0;
float GAP_EXTN_PENALTY = -0.1;

float align(int height, int width)
{
    register int i, j;
    float score;
    register float h, v, d;

    for (i = 0; i <= height; i++) {
        VAL[i][0] = 0;
        MOV[i][0] = 0;
    }

    for (j = 0; j <= width; j++) {
        VAL[0][j] = 0;
        MOV[0][j] = 0;
    }

    for (i = 1; i < height + 1; i++)
        for (j = 1; j < width + 1; j++) {
            h = VAL[i - 1][j];
            if (MOV[i - 1][j] != 1)
                h += GAP_OPEN_PENALTY;
            else
                h += GAP_EXTN_PENALTY;
            v = VAL[i][j - 1];
            if (MOV[i][j - 1] != 2)
                v += GAP_OPEN_PENALTY;
            else
                v += GAP_EXTN_PENALTY;
            d = VAL[i - 1][j - 1] + SCORE[i - 1][j - 1];
            if (d >= h && d >= v) {
                VAL[i][j] = d;
                MOV[i][j] = 0;
            }
            else if (h >= d && h >= v) {
                VAL[i][j] = h;
                MOV[i][j] = 1;
            }
            else {
                VAL[i][j] = v;
                MOV[i][j] = 2;
            }
        }

    score = VAL[height][width];

    return score;
}

int** alignment;

int process_alignment(int height, int width)
{
    int i, j, nali, id, ri, rj, n;
    char c;
    FILE* out;
    int ali = 0;

    i = height;
    j = width;
    ri = rj = 0;
    ali = 0;

    do {
        switch (MOV[i][j]) {
        case 0:
            i--;
            j--;
            break;
        case 1:
            i--;
            break;
        case 2:
            j--;
            break;
        }
        ali++;
    } while (i > 0 && j > 0);

    i = height;
    j = width;
    id = 0;
    ri = rj = 0;
    n = 0;

    do {
        switch (MOV[i][j]) {
        case 0:
            i--;
            j--;
            alignment[ali - n][0] = i;
            alignment[ali - n][1] = j;
            ri++;
            rj++;
            id++;
            break;
        case 1:
            i--;
            ri++;
            alignment[ali - n][0] = i;
            alignment[ali - n][1] = -1;
            break;
        case 2:
            j--;
            rj++;
            alignment[ali - n][0] = -1;
            alignment[ali - n][1] = j;
            break;
        }
        n++;
    } while (i > 0 && j > 0);

    return ali;
}

void vcross(float ax, float ay, float az, float bx, float by, float bz, float* cx, float* cy, float* cz)
{
    *cx = ay * bz - by * az;
    *cy = az * bx - bz * ax;
    *cz = ax * by - bx * ay;
}

float vdot(float ax, float ay, float az, float bx, float by, float bz)
{
    return ax * bx + ay * by + az * bz;
}

float calc_torsion(atom_struct* a1, atom_struct* a2, atom_struct* a3, atom_struct* a4)
{
    float v12x, v12y, v12z;
    float v43x, v43y, v43z;
    float zx, zy, zz;
    float px, py, pz;
    float xx, xy, xz;
    float yx, yy, yz;
    float u, v, angle;

    v12x = a1->x - a2->x;
    v12y = a1->y - a2->y;
    v12z = a1->z - a2->z;

    v43x = a4->x - a3->x;
    v43y = a4->y - a3->y;
    v43z = a4->z - a3->z;

    zx = a2->x - a3->x;
    zy = a2->y - a3->y;
    zz = a2->z - a3->z;

    vcross(zx, zy, zz, v12x, v12y, v12z, &px, &py, &pz);
    vcross(zx, zy, zz, v43x, v43y, v43z, &xx, &xy, &xz);
    vcross(zx, zy, zz, xx, xy, xz, &yx, &yy, &yz);

    u = vdot(xx, xy, xz, xx, xy, xz);
    v = vdot(yx, yy, yz, yx, yy, yz);

    angle = 360.;

    if (u < 0. || v < 0.)
        return angle;

    u = vdot(px, py, pz, xx, xy, xz) / sqrt(u);
    v = vdot(px, py, pz, yx, yy, yz) / sqrt(v);

    if (u != 0.0 || v != 0.0)
        angle = atan2(v, u) * (180.0f / M_PI);

    return angle;
}

int calc_chi(atom_struct* atom, float* chi, int which)
{
    atom_struct *a, *a1, *a2, *a3, *a4;
    int ri, atype;

    *chi = 0.0f;

    if (atom) {
        a = atom;
        atype = AA_NUMS[(int)a->aa];
        if (atype > 1) {
            ri = a->rnum;
            a1 = a2 = a3 = a4 = NULL;
            while (a && a->rnum == ri) {
                switch (which) {
                case 1:
                    if (nheavy[atype] > 1) {
                        if (a->aname[1] == 'N' && a->aname[2] == ' ')
                            a1 = a;
                        if (a->aname[1] == 'C' && a->aname[2] == 'A')
                            a2 = a;
                        if (a->aname[1] == heavy_atoms[10 * atype][0] && a->aname[2] == heavy_atoms[10 * atype][1] && a->aname[3] == heavy_atoms[10 * atype][2])
                            a3 = a;
                        if (a->aname[1] == heavy_atoms[10 * atype + 1][0] && a->aname[2] == heavy_atoms[10 * atype + 1][1] && a->aname[3] == heavy_atoms[10 * atype + 1][2])
                            a4 = a;
                    }
                    break;
                case 2:
                    if (nheavy[atype] > 2) {
                        if (a->aname[1] == 'C' && a->aname[2] == 'A')
                            a1 = a;
                        if (a->aname[1] == heavy_atoms[10 * atype][0] && a->aname[2] == heavy_atoms[10 * atype][1] && a->aname[3] == heavy_atoms[10 * atype][2])
                            a2 = a;
                        if (a->aname[1] == heavy_atoms[10 * atype + 1][0] && a->aname[2] == heavy_atoms[10 * atype + 1][1] && a->aname[3] == heavy_atoms[10 * atype + 1][2])
                            a3 = a;
                        if (a->aname[1] == heavy_atoms[10 * atype + 2][0] && a->aname[2] == heavy_atoms[10 * atype + 2][1] && a->aname[3] == heavy_atoms[10 * atype + 2][2])
                            a4 = a;
                    }
                    break;
                case 3:
                    if (nheavy[atype] > 3) {
                        if (a->aname[1] == heavy_atoms[10 * atype][0] && a->aname[2] == heavy_atoms[10 * atype][1] && a->aname[3] == heavy_atoms[10 * atype][2])
                            a1 = a;
                        if (a->aname[1] == heavy_atoms[10 * atype + 1][0] && a->aname[2] == heavy_atoms[10 * atype + 1][1] && a->aname[3] == heavy_atoms[10 * atype + 1][2])
                            a2 = a;
                        if (a->aname[1] == heavy_atoms[10 * atype + 2][0] && a->aname[2] == heavy_atoms[10 * atype + 2][1] && a->aname[3] == heavy_atoms[10 * atype + 2][2])
                            a3 = a;
                        if (a->aname[1] == heavy_atoms[10 * atype + 3][0] && a->aname[2] == heavy_atoms[10 * atype + 3][1] && a->aname[3] == heavy_atoms[10 * atype + 3][2])
                            a4 = a;
                    }
                    break;
                }
                a = a->next;
            }
            if (a1 && a2 && a3 && a4) {
                *chi = calc_torsion(a1, a2, a3, a4);
                return 1;
            }
        }
    }

    return 0;
}

int main(int argc, char** argv)
{
    FILE* out;
    mol_struct *mol1, *ligmol, *next_m;
    mol_struct* mol2;
    mol_struct *mol, *m1, *m2;
    atom_struct *atom, *a1, *a2, *lasta1, **alist1, **alist2, *ligatom;
    float **ca1, **ca2, **testca1, **testca2;
    float** stat;
    char *seq1, *seq2;
    int i, j, k, l, nca1, nca2, nali, m, n, o, all, first;
    float **rot, *trans, *cm;
    float rmsd, lrmsd;
    float x, y, z, dx, dy, dz, dd, tmpx, tmpy, tmpz;
    float cmx1, cmy1, cmz1;
    float cmx2, cmy2, cmz2;
    float cmr;
    char rid[6], rid2[6];
    float trmsd, srmsd;
    float dpos[2][100][3];
    float avg0, avg1, avg2;
    float sumsq0, sumsq1, sumsq2;
    float std0, std1, std2, dd1, dd2;
    float avg3, avg4, avg5;
    float sumsq3, sumsq4, sumsq5;
    float std3, std4, std5;
    float maxres;
    float z0, z1, z2, z3, z4, z5;
    float chi1_1, chi1_2, pchi1_1, pchi1_2, chid1;
    int ok1_1, ok1_2;
    float chi2_1, chi2_2, pchi2_1, pchi2_2, chid2;
    int ok2_1, ok2_2;
    float chi3_1, chi3_2, pchi3_1, pchi3_2, chid3;
    int ok3_1, ok3_2;

    printf("CCOMP Version 3.70 (c) 2002-2007 Piotr Rotkiewicz piotr@pirx.com\n");
    printf("Receptor/ligand complexes comparison program\n\n");

    setbuf(stdout, 0);

    for (i = 0; i < 255; i++)
        AA_NUMS[i] = 20;
    for (i = 0; i < 20; i++)
        AA_NUMS[(int)SHORT_AA_NAMES[i]] = i;

    if (argc < 2) {
        printf("Use: ccomp <file1.pdb> <file2.pdb> [all]\nBoth pdb1 and pdb2 files should include receptor (first) and ligand (second) molecules separated by TER line\n");
        printf("Use \"all\" as a third parameter for complete statistics (otherwise, only significant changes will be displayed (z-score>1)\n");
        exit(-1);
    }
    //printf("ok\n");

    printf("Input files: %s %s\n", argv[1], argv[2]);

    if (argc > 3 && !strcmp(argv[3], "all"))
        all = 1;
    else
        all = 0;

    mol1 = read_pdb(argv[1]);
    if (!mol1) {
        printf("ERROR: can't read the first molecule\n");
        exit(-2);
    }
    //printf("ok\n");

    mol2 = read_pdb(argv[2]);
    if (!mol2) {
        printf("ERROR: can't read the second molecule\n");
        exit(-3);
    }
    /*        
    if (mol1->natom!=mol2->natom) {
      printf("ERROR: number of atoms in both molecules doesn't match: %d %d\n", mol1->natom, mol2->natom);
      exit(-4);
    }
*/
    nca1 = 0;
    atom = mol1->atoms;
    while (atom) {
        if (!strncmp(atom->aname, " CA", 3))
            nca1++;
        atom = atom->next;
    }

    nca2 = 0;
    atom = mol2->atoms;
    while (atom) {
        if (!strncmp(atom->aname, " CA", 3))
            nca2++;
        atom = atom->next;
    }

    printf("Numbers of C-alpha atoms: %d %d\n", nca1, nca2);
    printf("Resolutions: %.2f %.2f\n", mol1->res, mol2->res);

    maxres = -1.0;
    if (mol1->res > maxres)
        maxres = mol1->res;
    if (mol2->res > maxres)
        maxres = mol2->res;

    if (!nca1 || !nca2) {
        printf("ERROR: no C-alpha atoms\n");
        exit(-5);
    }

    SCORE = (float**)malloc(sizeof(float*) * (nca1 + 1));
    for (i = 0; i < nca1 + 1; i++)
        SCORE[i] = (float*)malloc(sizeof(float) * (nca2 + 1));

    VAL = (float**)malloc(sizeof(float*) * (nca1 + 2));
    for (i = 0; i < nca1 + 2; i++)
        VAL[i] = (float*)malloc(sizeof(float) * (nca2 + 2));

    MOV = (char**)malloc(sizeof(char*) * (nca1 + 2));
    for (i = 0; i < nca1 + 2; i++)
        MOV[i] = (char*)malloc(sizeof(char) * (nca2 + 2));

    ca1 = (float**)calloc(sizeof(float*) * nca1, 1);
    for (i = 0; i < nca1; i++)
        ca1[i] = (float*)calloc(sizeof(float) * 3, 1);

    ca2 = (float**)calloc(sizeof(float*) * nca2, 1);
    for (i = 0; i < nca2; i++)
        ca2[i] = (float*)calloc(sizeof(float) * 3, 1);

    seq1 = (char*)calloc(sizeof(char) * (nca1 + 1), 1);
    seq2 = (char*)calloc(sizeof(char) * (nca2 + 1), 1);

    alist1 = (atom_struct**)calloc(sizeof(atom_struct*) * nca1, 1);
    alist2 = (atom_struct**)calloc(sizeof(atom_struct*) * nca2, 1);

    stat = (float**)calloc(sizeof(float*) * (nca1 + nca2), 1);
    for (i = 0; i < nca1 + nca2; i++)
        stat[i] = (float*)calloc(sizeof(float) * (10), 1);

    alignment = (int**)calloc(sizeof(int*) * (nca1 + nca2), 1);
    for (i = 0; i < nca1 + nca2; i++)
        alignment[i] = (int*)calloc(sizeof(int) * 2, 1);

    i = 0;
    atom = mol1->atoms;
    while (atom) {
        if (!strncmp(atom->aname, " CA", 3)) {
            ca1[i][0] = atom->x;
            ca1[i][1] = atom->y;
            ca1[i][2] = atom->z;
            seq1[i] = atom->aa;
            i++;
        }
        atom = atom->next;
    }

    i = 0;
    atom = mol2->atoms;
    while (atom) {
        if (!strncmp(atom->aname, " CA", 3)) {
            ca2[i][0] = atom->x;
            ca2[i][1] = atom->y;
            ca2[i][2] = atom->z;
            seq2[i] = atom->aa;
            i++;
        }
        atom = atom->next;
    }

    for (i = 0; i < nca1; i++)
        for (j = 0; j < nca2; j++)
            if (seq1[i] == seq2[j])
                SCORE[i][j] = 1.0;
            else
                SCORE[i][j] = 0.0;

    printf("\nAlignment score: %.1f\n", align(nca1, nca2));
    align(nca1, nca2);

    nali = process_alignment(nca1, nca2);
    printf("Alignment length: %d\n\n", nali);

    printf("Sequence 1:\n");

    j = 0;
    for (i = 0; i < nali; i++) {
        if (alignment[i][0] >= 0)
            printf("%c", seq1[alignment[i][0]]);
        else
            printf("-");
        if (j++ > 70) {
            printf("\n");
            j = 0;
        }
    }
    printf("\n\n");
    printf("Sequence 2:\n");
    j = 0;
    for (i = 0; i < nali; i++) {
        if (alignment[i][1] >= 0)
            printf("%c", seq2[alignment[i][1]]);
        else
            printf("-");
        if (j++ > 70) {
            printf("\n");
            j = 0;
        }
    }
    printf("\n\n");

    testca1 = (float**)calloc(sizeof(float*) * nali, 1);
    for (i = 0; i < nali; i++)
        testca1[i] = (float*)calloc(sizeof(float) * 3, 1);

    testca2 = (float**)calloc(sizeof(float*) * nali, 1);
    for (i = 0; i < nali; i++)
        testca2[i] = (float*)calloc(sizeof(float) * 3, 1);

    j = 0;
    for (i = 0; i < nali; i++) {
        if (alignment[i][0] >= 0 && alignment[i][1] >= 0 && seq1[alignment[i][0]] == seq2[alignment[i][1]]) {
            testca1[j][0] = ca1[alignment[i][0]][0];
            testca1[j][1] = ca1[alignment[i][0]][1];
            testca1[j][2] = ca1[alignment[i][0]][2];
            testca2[j][0] = ca2[alignment[i][1]][0];
            testca2[j][1] = ca2[alignment[i][1]][1];
            testca2[j][2] = ca2[alignment[i][1]][2];
            j++;
        }
    }

    printf("Length of aligned part: %d\n", j);

    rot = (float**)calloc(sizeof(float*) * 3, 1);
    for (i = 0; i < 3; i++)
        rot[i] = (float*)calloc(sizeof(float) * 3, 1);

    trans = (float*)calloc(sizeof(float) * 3, 1);
    cm = (float*)calloc(sizeof(float) * 3, 1);

    rmsd = superimpose(testca1, testca2, j, rot, trans, cm);

    printf("\nC-alpha  RMSD: %8.3f on %4d atoms\n\n", rmsd, j);

    // re-arrange molecules... tbd

    // transform all molecules

    mol = mol2;
    while (mol) {
        atom = mol->atoms;
        while (atom) {
            x = atom->x - cm[0];
            y = atom->y - cm[1];
            z = atom->z - cm[2];
            tmpx = x * rot[0][0] + y * rot[0][1] + z * rot[0][2];
            tmpy = x * rot[1][0] + y * rot[1][1] + z * rot[1][2];
            tmpz = x * rot[2][0] + y * rot[2][1] + z * rot[2][2];
            atom->x = tmpx + trans[0];
            atom->y = tmpy + trans[1];
            atom->z = tmpz + trans[2];
            atom = atom->next;
        }
        mol = mol->next;
    }

    // compare molecules

    m1 = mol1;
    m2 = mol2;
    j = 0;
    rid[5] = 0;
    while (m1 && m2) {
        trmsd = srmsd = 0.0;
        lrmsd = 0.0;
        i = 0;
        k = 0;
        a1 = m1->atoms;
        a2 = m2->atoms;
        strncpy(rid, a1->rid, 5);
        lasta1 = NULL;

        if (m1 == mol1 && m2 == mol2) { // first pair: receptors
            a1 = m1->atoms;
            i = 0;
            strncpy(rid, "xxxxx", 5);
            do {
                if (strncmp(rid, a1->rid, 5)) {
                    alist1[i] = a1;
                    //printf("a1list %d %d\n", i, a1);
                    i++;
                    strncpy(rid, a1->rid, 5);
                }
                a1 = a1->next;
            } while (a1);

            a2 = m2->atoms;
            i = 0;
            strncpy(rid, "xxxxx", 5);
            do {
                if (strncmp(rid, a2->rid, 5)) {
                    alist2[i] = a2;
                    //printf("a2list %d %d\n", i, a2);
                    i++;
                    strncpy(rid, a2->rid, 5);
                }
                a2 = a2->next;
            } while (a2);
            //printf("a1,a2\n");
            l = 0;

            printf("Molecule pair %2d:\n\n", j + 1);

            for (i = 0; i < nali; i++) {
                if (alignment[i][0] >= 0 && alignment[i][1] >= 0 && seq1[alignment[i][0]] == seq2[alignment[i][1]]) {
                    a1 = alist1[alignment[i][0]];
                    a2 = alist2[alignment[i][1]];

                    strncpy(rid, a1->rid, 5);
                    strncpy(rid2, a2->rid, 5);

                    lrmsd = 0.0;

                    k = 0;

                    cmx1 = cmy1 = cmz1 = 0.0;
                    cmx2 = cmy2 = cmz2 = 0.0;

                    do {
                        if (a1->aname[0] != 'H' && a2->aname[0] != 'H' && !(a1->aname[0] == ' ' && a1->aname[1] == 'H') && !(a2->aname[0] == ' ' && a2->aname[1] == 'H')) {
                            dpos[0][k][0] = a1->x;
                            dpos[0][k][1] = a1->y;
                            dpos[0][k][2] = a1->z;
                            dpos[1][k][0] = a2->x;
                            dpos[1][k][1] = a2->y;
                            dpos[1][k][2] = a2->z;
                            dx = a1->x - a2->x;
                            dy = a1->y - a2->y;
                            dz = a1->z - a2->z;
                            cmx1 += a1->x;
                            cmy1 += a1->y;
                            cmz1 += a1->z;
                            cmx2 += a2->x;
                            cmy2 += a2->y;
                            cmz2 += a2->z;
                            trmsd += dx * dx + dy * dy + dz * dz;
                            lrmsd += dx * dx + dy * dy + dz * dz;
                            k++;
                            l++;
                        }
                        a1 = a1->next;
                        a2 = a2->next;
                    } while (a1 && a2 && !strncmp(rid, a1->rid, 5) && !strncmp(rid2, a2->rid, 5));

                    if (k) {
                        cmx1 /= (float)k;
                        cmy1 /= (float)k;
                        cmz1 /= (float)k;
                        cmx2 /= (float)k;
                        cmy2 /= (float)k;
                        cmz2 /= (float)k;
                    }

                    stat[i][0] = k;
                    if (k)
                        stat[i][1] = sqrt(lrmsd / (float)k);
                    else
                        stat[i][1] = 0.0; // LOCAL RMSD
                    stat[i][2] = sqrt((cmx1 - cmx2) * (cmx1 - cmx2) + (cmy1 - cmy2) * (cmy1 - cmy2) + (cmz1 - cmz2) * (cmz1 - cmz2));

                    o = 0;
                    lrmsd = 0.0;

                    for (m = 0; m < k - 1; m++) {
                        for (n = m + 1; n < k; n++) {
                            cmx1 = dpos[0][m][0] - dpos[0][n][0];
                            cmx2 = dpos[1][m][0] - dpos[1][n][0];
                            cmy1 = dpos[0][m][1] - dpos[0][n][1];
                            cmy2 = dpos[1][m][1] - dpos[1][n][1];
                            cmz1 = dpos[0][m][2] - dpos[0][n][2];
                            cmz2 = dpos[1][m][2] - dpos[1][n][2];
                            dd1 = sqrt(cmx1 * cmx1 + cmy1 * cmy1 + cmz1 * cmz1);
                            dd2 = sqrt(cmx2 * cmx2 + cmy2 * cmy2 + cmz2 * cmz2);
                            lrmsd += (dd1 - dd2) * (dd1 - dd2);
                            o++;
                        }
                    }

                    if (o) {
                        stat[i][3] = sqrt(lrmsd / (float)o); // DISTANCE RMSD
                    }
                    else {
                        stat[i][3] = 0.0;
                    }

                    lrmsd = 0.0;
                    k = 0;
                    cmx1 = cmy1 = cmz1 = 0.0;
                    cmx2 = cmy2 = cmz2 = 0.0;

                    a1 = alist1[alignment[i][0]];
                    a2 = alist2[alignment[i][1]];
                    do {
                        if (a1->flags & FLAG_SIDECHAIN && a2->flags & FLAG_SIDECHAIN) {
                            if (a1->aname[0] != 'H' && a2->aname[0] != 'H' && !(a1->aname[0] == ' ' && a1->aname[1] == 'H') && !(a2->aname[0] == ' ' && a2->aname[1] == 'H')) {
                                //printf("%d %s\n", i, a1->aname);
                                dpos[0][k][0] = a1->x;
                                dpos[0][k][1] = a1->y;
                                dpos[0][k][2] = a1->z;
                                dpos[1][k][0] = a2->x;
                                dpos[1][k][1] = a2->y;
                                dpos[1][k][2] = a2->z;
                                dx = a1->x - a2->x;
                                dy = a1->y - a2->y;
                                dz = a1->z - a2->z;
                                cmx1 += a1->x;
                                cmy1 += a1->y;
                                cmz1 += a1->z;
                                cmx2 += a2->x;
                                cmy2 += a2->y;
                                cmz2 += a2->z;
                                srmsd += dx * dx + dy * dy + dz * dz;
                                lrmsd += dx * dx + dy * dy + dz * dz;
                                k++;
                            }
                        }
                        a1 = a1->next;
                        a2 = a2->next;
                    } while (a1 && a2 && !strncmp(rid, a1->rid, 5) && !strncmp(rid2, a2->rid, 5));

                    if (k) {
                        cmx1 /= (float)k;
                        cmy1 /= (float)k;
                        cmz1 /= (float)k;

                        cmx2 /= (float)k;
                        cmy2 /= (float)k;
                        cmz2 /= (float)k;

                        stat[i][4] = sqrt(lrmsd / (float)k); // LOCAL RMSD
                        stat[i][5] = sqrt((cmx1 - cmx2) * (cmx1 - cmx2) + (cmy1 - cmy2) * (cmy1 - cmy2) + (cmz1 - cmz2) * (cmz1 - cmz2));
                    }

                    o = 0;
                    lrmsd = 0.0;
                    if (k > 1) {
                        for (m = 0; m < k - 1; m++) {
                            for (n = m + 1; n < k; n++) {
                                cmx1 = dpos[0][m][0] - dpos[0][n][0];
                                cmx2 = dpos[1][m][0] - dpos[1][n][0];
                                cmy1 = dpos[0][m][1] - dpos[0][n][1];
                                cmy2 = dpos[1][m][1] - dpos[1][n][1];
                                cmz1 = dpos[0][m][2] - dpos[0][n][2];
                                cmz2 = dpos[1][m][2] - dpos[1][n][2];
                                dd1 = sqrt(cmx1 * cmx1 + cmy1 * cmy1 + cmz1 * cmz1);
                                dd2 = sqrt(cmx2 * cmx2 + cmy2 * cmy2 + cmz2 * cmz2);
                                lrmsd += (dd1 - dd2) * (dd1 - dd2);
                                o++;
                            }
                        }
                        if (o) {
                            stat[i][6] = sqrt(lrmsd / (float)o); // DISTANCE RMSD
                        }
                        else {
                            stat[i][6] = 0.0;
                        }
                    }
                    else {
                        stat[i][6] = 0.0; // DISTANCE RMSD
                    }
                }
            }

            printf("Calculation completed for the receptor; %d amino acids aligned\n\n", nali);

            avg0 = avg1 = avg2 = avg3 = avg4 = avg5 = 0.0;
            sumsq0 = sumsq1 = sumsq2 = sumsq3 = sumsq4 = sumsq5 = 0.0;

            o = 0;
            for (i = 0; i < nali; i++) {
                if (alignment[i][0] >= 0 && alignment[i][1] >= 0 && seq1[alignment[i][0]] == seq2[alignment[i][1]]) {
                    avg0 += stat[i][1];
                    avg1 += stat[i][2];
                    avg2 += stat[i][3];
                    avg3 += stat[i][4];
                    avg4 += stat[i][5];
                    avg5 += stat[i][6];
                    sumsq0 += stat[i][1] * stat[i][1];
                    sumsq1 += stat[i][2] * stat[i][2];
                    sumsq2 += stat[i][3] * stat[i][3];
                    sumsq3 += stat[i][4] * stat[i][4];
                    sumsq4 += stat[i][5] * stat[i][5];
                    sumsq5 += stat[i][6] * stat[i][6];
                    o++;
                }
            }

            if (o) {
                avg0 /= (float)o;
                avg1 /= (float)o;
                avg2 /= (float)o;
                avg3 /= (float)o;
                avg4 /= (float)o;
                avg5 /= (float)o;
                sumsq0 /= (float)o;
                sumsq1 /= (float)o;
                sumsq2 /= (float)o;
                sumsq3 /= (float)o;
                sumsq4 /= (float)o;
                sumsq5 /= (float)o;
            }

            std0 = sumsq0 - avg0 * avg0;
            std1 = sumsq1 - avg1 * avg1;
            std2 = sumsq2 - avg2 * avg2;
            std3 = sumsq3 - avg3 * avg3;
            std4 = sumsq4 - avg4 * avg4;
            std5 = sumsq5 - avg5 * avg5;

            printf("C-alpha atoms are superimposed according to the sequence alignment.\n");
            printf("L-RMSD: local all-atom RMSD per residuum: reflects both spatial differences and conformational changes\n");
            printf("C-RMSD: center of mass deviation per residuum: reflects mostly spatial differences\n");
            printf("D-RMSD: distance RMSD per residuum: reflects only conformational changes\n\n");
            printf("-Z: z-score (normalized deviation, z-score = score/sigma)\n");
            printf("+*#: mark z-score>1.0, for L, C, D RMSD, respectively\n");

            printf("\nAverage  RMSD (whole aa): %10.5f, CM-RMSD: %10.5f, DISTANCE-RMSD: %10.5f\n", avg0, avg1, avg2);
            printf("Std.dev. RMSD (whole aa): %10.5f, CM-RMSD: %10.5f, DISTANCE-RMSD: %10.5f\n", std0, std1, std2);
            printf("\nAverage  RMSD (side chain): %10.5f, CM-RMSD: %10.5f, DISTANCE-RMSD: %10.5f\n", avg3, avg4, avg5);
            printf("Std.dev. RMSD (side chain): %10.5f, CM-RMSD: %10.5f, DISTANCE-RMSD: %10.5f\n\n", std3, std4, std5);

            printf("                                        -------------------- WHOLE AMINO ACID ------------------        --------------------- SIDE CHAIN ONLY ------------------\n");
            printf("                                        L-RMSD  L-RMSD-Z    C-RMSD  C-RMSD-Z    D-RMSD  D-RMSD-Z        L-RMSD  L-RMSD-Z    C-RMSD  C-RMSD-Z    D-RMSD  D-RMSD-Z\n");

            first = 1;

            if (std2 < 1e-15)
                std2 = 1.0f;
            if (std5 < 1e-15)
                std5 = 1.0f;

            for (i = 0; i < nali; i++) {
                if (alignment[i][0] >= 0 && alignment[i][1] >= 0 && seq1[alignment[i][0]] == seq2[alignment[i][1]]) {
                    if (first) {
                        a1 = alist1[alignment[i][0]];
                        a2 = alist2[alignment[i][1]];
                        z0 = -(avg0 - stat[i][1]) / std0;
                        z1 = -(avg1 - stat[i][2]) / std1;
                        z2 = -(avg2 - stat[i][3]) / std2;
                        z3 = -(avg3 - stat[i][4]) / std3;
                        z4 = -(avg4 - stat[i][5]) / std4;
                        z5 = -(avg5 - stat[i][6]) / std5;
                        if (z0 > 1.0f || z1 > 1.0f || z2 > 1.0f || z3 > 1.0f || z4 > 1.0f || z5 > 1.0f || all) {
                            printf("%3s %5s - %3s %5s", a1->rname, a1->rid, a2->rname, a2->rid);
                            strncpy(rid, a1->rid, 5);
                            printf(": %4d atoms :   %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  ", (int)stat[i][0], stat[i][1], -(avg0 - stat[i][1]) / std0, stat[i][2], -(avg1 - stat[i][2]) / std1, stat[i][3], -(avg2 - stat[i][3]) / std2);
                            if (z0 > 1.0)
                                printf("+");
                            else
                                printf(" ");
                            if (z1 > 1.0)
                                printf("*");
                            else
                                printf(" ");
                            if (z2 > 1.0)
                                printf("#");
                            else
                                printf(" ");
                            printf(" %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  ", stat[i][4], -(avg3 - stat[i][4]) / std3, stat[i][5], -(avg4 - stat[i][5]) / std4, stat[i][6], -(avg5 - stat[i][6]) / std5);
                            if (z3 > 1.0)
                                printf("+");
                            else
                                printf(" ");
                            if (z4 > 1.0)
                                printf("*");
                            else
                                printf(" ");
                            if (z5 > 1.0)
                                printf("#");
                            else
                                printf(" ");
                            printf("\n");
                        }
                    }
                    first = 1;
                }
            }
            //        printf("\nSummary: %4d atoms   Average RMSD per residuum: %8.3f (whole),  %8.3f (side chain)\n", l, sqrt(trmsd/(float)l), sqrt(srmsd/(float)l));

            printf("\n\nChi-angles analysis");
            if (!all)
                printf(" (only chi differences > 10 deg are displayed)\n");
            else
                printf("\n");

            printf("\n");
            printf("                           CHI(1)  CHI(1)   CHI(1)       CHI(2)  CHI(2)   CHI(2)       CHI(3)  CHI(3)   CHI(3)\n");
            printf("                           FIRST   SECOND  DIFFER.       FIRST   SECOND  DIFFER.       FIRST   SECOND  DIFFER.\n");
            printf("                         -------  -------  -------     -------  -------  -------     -------  -------  -------\n");

            for (i = 0; i < nali; i++) {
                if (alignment[i][0] >= 0 && alignment[i][1] >= 0 && seq1[alignment[i][0]] == seq2[alignment[i][1]]) {
                    if (first) {
                        a1 = alist1[alignment[i][0]];
                        a2 = alist2[alignment[i][1]];

                        ok1_1 = calc_chi(a1, &chi1_1, 1);
                        ok1_2 = calc_chi(a2, &chi1_2, 1);
                        pchi1_1 = chi1_1;
                        if (pchi1_1 < 0.0)
                            pchi1_1 += 360.0f;
                        pchi1_2 = chi1_2;
                        if (pchi1_2 < 0.0)
                            pchi1_2 += 360.0f;
                        chid1 = pchi1_1 - pchi1_2;
                        if (chid1 > 180.0f)
                            chid1 = 360.0f - chid1;
                        if (chid1 < -180.0f)
                            chid1 = 360.0f + chid1;

                        ok2_1 = calc_chi(a1, &chi2_1, 2);
                        ok2_2 = calc_chi(a2, &chi2_2, 2);
                        pchi2_1 = chi2_1;
                        if (pchi2_1 < 0.0)
                            pchi2_1 += 360.0f;
                        pchi2_2 = chi2_2;
                        if (pchi2_2 < 0.0)
                            pchi2_2 += 360.0f;
                        chid2 = pchi2_1 - pchi2_2;
                        if (chid2 > 180.0f)
                            chid2 = 360.0f - chid2;
                        if (chid2 < -180.0f)
                            chid2 = 360.0f + chid2;

                        ok3_1 = calc_chi(a1, &chi3_1, 3);
                        ok3_2 = calc_chi(a2, &chi3_2, 3);
                        pchi3_1 = chi3_1;
                        if (pchi3_1 < 0.0)
                            pchi3_1 += 360.0f;
                        pchi3_2 = chi3_2;
                        if (pchi3_2 < 0.0)
                            pchi3_2 += 360.0f;
                        chid3 = pchi3_1 - pchi3_2;
                        if (chid3 > 180.0f)
                            chid3 = 360.0f - chid3;
                        if (chid3 < -180.0f)
                            chid3 = 360.0f + chid3;

                        if (all || (ok1_1 && ok1_2 && fabs(chid1) > 10.0f) || (ok2_1 && ok2_2 && fabs(chid2) > 10.0f) || (ok3_1 && ok3_2 && fabs(chid3) > 10.0f)) {
                            printf("%3s %5s - %3s %5s : ", a1->rname, a1->rid, a2->rname, a2->rid);
                            if (ok1_1 && ok1_2) {
                                printf("%8.3f %8.3f %8.3f    ", chi1_1, chi1_2, chid1);
                            }
                            else {
                                printf("%8s %8s %8s    ", "-", "-", "-");
                            }
                            if (ok2_1 && ok2_2) {
                                printf("%8.3f %8.3f %8.3f    ", chi2_1, chi2_2, chid2);
                            }
                            else {
                                printf("%8s %8s %8s    ", "-", "-", "-");
                            }
                            if (ok3_1 && ok3_2) {
                                printf("%8.3f %8.3f %8.3f    ", chi3_1, chi3_2, chid3);
                            }
                            else {
                                printf("%8s %8s %8s    ", "-", "-", "-");
                            }
                            printf("\n");
                        }
                    }
                }
            }
        }
        else {
            printf("\nMolecule pair %2d:\n\n", j + 1);
            i = k = 0;
            rmsd = lrmsd = 0.0;
            cmx1 = cmy1 = cmz1 = cmx2 = cmy2 = cmz2 = 0.0f;
            while (a1 && a2) {
                if (k && lasta1 && strncmp(rid, a1->rid, 5)) { // next residuum
                    printf("%3s %5s : %4d atoms : local RMSD %8.3f, CM-RMSD: %8.3f\n", lasta1->rname, lasta1->rid, k, sqrt(lrmsd / (float)k), sqrt((cmx1 - cmx2) * (cmx1 - cmx2) + (cmy1 - cmy2) * (cmy1 - cmy2) + (cmz1 - cmz2) * (cmz1 - cmz2)));
                    lrmsd = 0.0;
                    k = 0;
                    strncpy(rid, a1->rid, 5);
                }
                dx = a1->x - a2->x;
                dy = a1->y - a2->y;
                dz = a1->z - a2->z;
                cmx1 += a1->x;
                cmy1 += a1->y;
                cmz1 += a1->z;
                cmx2 += a2->x;
                cmy2 += a2->y;
                cmz2 += a2->z;
                rmsd += dx * dx + dy * dy + dz * dz;
                lrmsd += dx * dx + dy * dy + dz * dz;
                lasta1 = a1;
                a1 = a1->next;
                a2 = a2->next;
                i++;
                k++;
            }
            if (k && lasta1) { // last residuum
                cmx1 /= k;
                cmy1 /= k;
                cmz1 /= k;
                cmx2 /= k;
                cmy2 /= k;
                cmz2 /= k;
                printf("%3s %5s : %4d atoms : local RMSD %8.3f, CM-RMSD: %8.3f\n", lasta1->rname, lasta1->rid, k, sqrt(lrmsd / (float)k), sqrt((cmx1 - cmx2) * (cmx1 - cmx2) + (cmy1 - cmy2) * (cmy1 - cmy2) + (cmz1 - cmz2) * (cmz1 - cmz2)));
                lrmsd = 0.0;
                k = 0;
            }
            //        printf("\nSummary: %4d atoms   Average RMSD per residuum: %8.3f (whole),  %8.3f (side chain)\n", l, sqrt(trmsd/(float)l), sqrt(srmsd/(float)l));
        }
        m1 = m1->next;
        m2 = m2->next;
        j++;
    }

    printf("\nWriting output file: complex_superimposed.pdb\n");

    out = fopen("complex_superimposed.pdb", "w");

    mol = mol1;
    while (mol) {
        atom = mol->atoms;
        while (atom) {
            if (atom->het)
                fprintf(out, "HETATM");
            else
                fprintf(out, "ATOM  ");
            fprintf(out, "%5d", atom->anum);
            fprintf(out, " %4s ", atom->aname);
            fprintf(out, "%3s ", atom->rname);
            fprintf(out, "%c", atom->chain);
            fprintf(out, "%5s   ", atom->rid);
            fprintf(out, "%8.3f%8.3f%8.3f", atom->x, atom->y, atom->z);
            fprintf(out, "%6.2f%6.2f\n", atom->dat1, atom->dat2);
            atom = atom->next;
        }
        fprintf(out, "TER\n");
        mol = mol->next;
    }

    mol = mol2;
    while (mol) {
        atom = mol->atoms;
        while (atom) {
            if (atom->het)
                fprintf(out, "HETATM");
            else
                fprintf(out, "ATOM  ");
            fprintf(out, "%5d", atom->anum);
            fprintf(out, " %4s ", atom->aname);
            fprintf(out, "%3s ", atom->rname);
            fprintf(out, "%c", atom->chain);
            fprintf(out, "%5s   ", atom->rid);
            fprintf(out, "%8.3f%8.3f%8.3f", atom->x, atom->y, atom->z);
            fprintf(out, "%6.2f%6.2f\n", atom->dat1, atom->dat2);
            atom = atom->next;
        }
        fprintf(out, "TER\n");
        mol = mol->next;
    }
    fprintf(out, "END\n");

    fclose(out);

    printf("\nReceptor-ligand contact analysis\n");
    printf("Displays all atoms < 3.5A\n");

    printf("\n--- %s ---\n", argv[1]);

    mol = mol1;
    if (mol && mol->next) {
        ligmol = mol->next;
        atom = mol->atoms;
        sprintf(rid, "xxx");
        while (atom) {
            ligatom = ligmol->atoms;
            while (ligatom) {
                dx = atom->x - ligatom->x;
                dy = atom->y - ligatom->y;
                dz = atom->z - ligatom->z;
                dd = dx * dx + dy * dy + dz * dz;
                if (dd < 3.5 * 3.5) {
                    if (strncmp(rid, atom->rid, 3)) {
                        strncpy(rid, atom->rid, 3);
                        printf("\n%3s%3s\n", atom->rname, atom->rid);
                    }
                    printf("   %s[%4d]  ---  %s %s[%4d] : %6.2f A\n", atom->aname, atom->anum, ligatom->rname, ligatom->aname, ligatom->anum, sqrt(dd));
                }
                ligatom = ligatom->next;
            }
            atom = atom->next;
        }
    }

    printf("\n\n--- %s ---\n", argv[2]);

    mol = mol2;
    if (mol && mol->next) {
        ligmol = mol->next;
        atom = mol->atoms;
        sprintf(rid, "xxx");
        while (atom) {
            ligatom = ligmol->atoms;
            while (ligatom) {
                dx = atom->x - ligatom->x;
                dy = atom->y - ligatom->y;
                dz = atom->z - ligatom->z;
                dd = dx * dx + dy * dy + dz * dz;
                if (dd < 3.5 * 3.5) {
                    if (strncmp(rid, atom->rid, 3)) {
                        strncpy(rid, atom->rid, 3);
                        printf("\n%3s%3s\n", atom->rname, atom->rid);
                    }
                    printf("   %s[%4d]  ---  %s %s[%4d] : %6.2f A\n", atom->aname, atom->anum, ligatom->rname, ligatom->aname, ligatom->anum, sqrt(dd));
                }
                ligatom = ligatom->next;
            }
            atom = atom->next;
        }
    }

    printf("\nCCOMP done.\n");

    for (i = 0; i < nca1 + 1; i++)
        free(SCORE[i]);
    free(SCORE);

    atom = mol1->atoms;
    while (atom) {
        a1 = atom->next;
        free(atom);
        atom = a1;
    }

    while (mol1) {
        next_m = mol1->next;
        free(mol1);
        mol1 = next_m;
    }

    atom = mol2->atoms;
    while (atom) {
        a1 = atom->next;
        free(atom);
        atom = a1;
    }

    while (mol2) {
        next_m = mol2->next;
        free(mol2);
        mol2 = next_m;
    }

    free(alist1);
    free(alist2);

    return 0;
}
