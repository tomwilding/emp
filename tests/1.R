test.examples <- function()
{
  checkEquals(getObservations(c(1,2,3,4,5), 0.5),c(1,3,5))
}
 
test.deactivation <- function()
{
  DEACTIVATED('Deactivating this test function')
}