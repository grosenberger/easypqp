import unittest
from easypqp.convert import get_scan


class TestConvert(unittest.TestCase):

    def test_get_scan(self):
        self.assertEqual(11, get_scan("controllerType=0 controllerNumber=1 scan=2 demux=0", 11))
        self.assertEqual(11, get_scan("sample=2 period=3 cycle=4 experiment=5", 11))
        self.assertEqual(11, get_scan("frame=2 scan=3", 11))

        self.assertEqual(11, get_scan("controllerType=0 controllerNumber=1 scan=11", 22))
        self.assertEqual(11, get_scan("function=0 process=1 scan=11", 22))
        self.assertEqual(11, get_scan("jobRun=0 spotLabel=asw spectrum=11", 22))
        self.assertEqual(11, get_scan("11", 22))
        self.assertEqual(11, get_scan("scan=11", 22))
        self.assertEqual(11, get_scan("spectrum=11", 22))
        self.assertEqual(11, get_scan("scanId=11", 22))
        self.assertEqual(11, get_scan("index=11", 22))
        self.assertEqual(11, get_scan("frame=11", 22))


if __name__ == '__main__':
    unittest.main()



