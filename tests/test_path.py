from unittest import TestCase

import jade.basic.path as jade_path

class TestPath(TestCase):
    def test_find_database(self):
        s = jade_path.get_database_path()
        self.assertTrue(isinstance(s, str))

        contents = jade_path.parse_contents(jade_path.get_database_testing_path())
        self.assertEquals(contents[0], "TEST_ASSERT")