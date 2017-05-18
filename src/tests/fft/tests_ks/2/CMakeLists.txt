ADDMPITEST(ks test_ks_2 2)

add_test(check_ks_2 python check_output.py)
set_tests_properties(check_ks_2 PROPERTIES DEPENDS test_ks_2)
