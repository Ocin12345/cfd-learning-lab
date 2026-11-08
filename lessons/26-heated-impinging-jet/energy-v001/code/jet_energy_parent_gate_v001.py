"""Repeat parent/source/resource admission after waits and final races."""
def wait_checked(check_parent,gate,sleep):
 while True:
  check_parent()
  try:
   gate();check_parent();return gate()
  except BlockingIOError:sleep(15)
