
#include <stdio.h>
#include <stdlib.h>
#include "libaxb.h"


#define AXB_ERR_CHECK( arg )   do { int err_check_12388123 = arg; if (err_check_12388123) { fprintf(stderr, "Failure in line %d in file %s\n", __LINE__, __FILE__); return err_check_12388123; }  } while (0)


int main(int argc, char **argv)
{
  axbHandle_t axbHandle;

  printf("Initializing libaxb...\n");
  AXB_ERR_CHECK(axbInit(&axbHandle));

  double ones[] = {1, 1, 1, 1, 1};
  axbVec_t x;
  AXB_ERR_CHECK(axbVecCreateBegin(axbHandle, &x));
  AXB_ERR_CHECK(axbVecSetSize(x, 5));
  AXB_ERR_CHECK(axbVecCreateEnd(x));
  AXB_ERR_CHECK(axbVecSetValues(x, ones, AXB_REAL_DOUBLE));

  double twos[] = {2, 2, 2, 2, 2};
  axbVec_t y;
  AXB_ERR_CHECK(axbVecCreateBegin(axbHandle, &y));
  AXB_ERR_CHECK(axbVecSetSize(y, 5));
  AXB_ERR_CHECK(axbVecCreateEnd(y));
  AXB_ERR_CHECK(axbVecSetValues(y, twos, AXB_REAL_DOUBLE));

  double three = 3.0;
  axbScalar_t alpha;
  AXB_ERR_CHECK(axbScalarCreate(axbHandle, &alpha, &three, AXB_REAL_DOUBLE, NULL));
  AXB_ERR_CHECK(axbVecAXPY(y, alpha, x));

  double result[] = {0, 0, 0, 0, 0};
  printf("Result before: %g %g %g %g %g\n", result[0], result[1], result[2], result[3], result[4]);
  AXB_ERR_CHECK(axbVecGetValues(y, result, AXB_REAL_DOUBLE));
  printf("Result after: %g %g %g %g %g\n", result[0], result[1], result[2], result[3], result[4]);

  printf("Shutting down libaxb...\n");
  AXB_ERR_CHECK(axbVecDestroy(x));
  AXB_ERR_CHECK(axbVecDestroy(y));
  AXB_ERR_CHECK(axbScalarDestroy(alpha));
  AXB_ERR_CHECK(axbFinalize(axbHandle));

  printf("Completed successfully\n");

  return EXIT_SUCCESS;
}
