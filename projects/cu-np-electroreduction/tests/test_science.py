"""SYNTHETIC SOFTWARE TESTS. These numbers are not scientific observations."""
import math
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from science import *

class ScientificChecks(unittest.TestCase):
    def test_unit_constants(self):
        self.assertAlmostEqual(HARTREE_EV*EV_KJMOL,2625.49964,places=4)
    def test_parse_complete_not_partial(self):
        self.assertEqual(parse_xtb("TOTAL ENERGY -12.25 Eh\nnormal termination of xtb")["energy_hartree"],-12.25)
        for text in ["TOTAL ENERGY -12.25 Eh", "TOTAL ENERGY -12.25 Eh\nnormal termination of xtb\nSCC is not converged"]:
            with self.assertRaises(ValueError): parse_xtb(text)
        with self.assertRaises(ValueError): parse_xtb("TOTAL ENERGY -12.25 Eh\nnormal termination of xtb","opt")
    def test_balance_electron(self):
        species={"H+":{"H":1,"charge":1},"e-":{"charge":-1},"H2":{"H":2,"charge":0}}
        self.assertTrue(balanced({"H+":-2,"e-":-2,"H2":1},species)[0])
        self.assertFalse(balanced({"H+":-2,"e-":-1,"H2":1},species)[0])
    def test_reference_offset_sign_and_guard(self):
        a=dict(solvent="test",temperature_K=298.15,common_reference="X",offset_V=.2,source="synthetic fixture")
        b={**a,"offset_V":.5}
        self.assertAlmostEqual(convert_potential(1.,a,b),.7)
        with self.assertRaises(ValueError): convert_potential(1.,a,{**b,"solvent":"water"})
        with self.assertRaises(ValueError): convert_potential(1.,{},b)
    def test_reservoirs(self):
        self.assertAlmostEqual(reservoir_free_energy(1.,{"H":1},{"H":.3}),.7)
        with self.assertRaises(ValueError): reservoir_free_energy(1.,{"H":1},{})
    def test_detailed_balance(self):
        f,r=rate_pair(0.,.1,.6)
        self.assertAlmostEqual(math.log(f/r),-.1/(KB_EV*298.15),places=11)
        with self.assertRaises(ValueError): rate_pair(0,.3,.2)
    def test_site_equilibrium_and_driving(self):
        theta,j=site_steady_state([1,1,1,1,1,1])
        self.assertAlmostEqual(sum(theta),1.)
        self.assertAlmostEqual(j,0.)
        _,j=site_steady_state([2,1,2,1,2,1])
        self.assertGreater(j,0)
    def test_energy_and_fe_recovery(self):
        self.assertAlmostEqual(electrical_kwh_kg(3,2,.8,200,.9),1.11672838106)
        with self.assertRaises(ValueError): electrical_kwh_kg(3,2,80,200)
    def test_no_split_leakage(self):
        rows=[{"identity":"a","group":"A"},{"identity":"a","group":"A"},{"identity":"b","group":"B"}]
        train,hold=grouped_split(rows,{"B"})
        self.assertEqual(len(train),2); self.assertEqual(len(hold),1)
        with self.assertRaises(ValueError): grouped_split(rows+[{"identity":"a","group":"B"}],{"B"})

if __name__ == "__main__": unittest.main()
