ADDMPITEST(ks test_ks_1 2)

add_test(check_ks_1 python check_output.py)
set_tests_properties(check_ks_1 PROPERTIES DEPENDS test_ks_1)
