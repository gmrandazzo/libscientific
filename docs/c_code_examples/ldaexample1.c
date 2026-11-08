/* Linear Discriminant Analysis: classify objects of 3 classes labelled 1, 2, 3. */
#include <stdio.h>
#include <scientific.h>

int main(void)
{
    matrix *x, *y, *features, *scores, *mnpdf, *labels;
    LDAMODEL *model;
    size_t i;

    NewMatrix(&x, 12, 2);
    NewMatrix(&y, 12, 1);
    for(i = 0; i < 12; i++){
        size_t c = i / 4;  /* 4 objects per class */
        x->data[i][0] = 10.0 * (double)c + 0.1 * (double)(i % 4);
        x->data[i][1] = -5.0 * (double)c + 0.2 * (double)((i * 3) % 4);
        y->data[i][0] = (double)(c + 1);
    }

    NewLDAModel(&model);
    LDA(x, y, model);
    printf("%zu classes, first label %zu\n", model->nclass, model->class_start);

    initMatrix(&features);
    initMatrix(&scores);
    initMatrix(&mnpdf);
    initMatrix(&labels);
    LDAPrediction(x, model, features, scores, mnpdf, labels);
    for(i = 0; i < x->row; i += 4)
        printf("object %zu: true class %g, predicted %g\n", i, y->data[i][0], labels->data[i][0]);

    DelMatrix(&labels);
    DelMatrix(&mnpdf);
    DelMatrix(&scores);
    DelMatrix(&features);
    DelLDAModel(&model);
    DelMatrix(&y);
    DelMatrix(&x);
    return 0;
}
