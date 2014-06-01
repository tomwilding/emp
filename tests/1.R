test.examples <- function()
{
  checkEquals(getObservations(c(1,2,3,4,5), 0.5), c(1,3,5))

  checkEquals(getObservations(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21), 0.05), c(1,21))

  checkEquals(reduceParams(c(1,1,1,1,1,1,1,1), 4), c(1,1,1,1))
}
 
test.deactivation <- function()
{
  DEACTIVATED('Deactivating this test function')
}