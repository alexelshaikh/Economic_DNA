import unittest

import storage_service


class ReferenceModelTests(unittest.TestCase):
    def test_zero_duration_has_only_initial_write(self):
        model = storage_service.amazon(n=1, k=1, obj_size_mb=100, start_year=2025, d=0)
        self.assertAlmostEqual(float(model.eval().doit()), float(model._get_initial_write_cost()))

    def test_replacement_at_end_of_horizon_is_excluded(self):
        model = storage_service.tape_on_premise(
            n=1, k=0, obj_size_mb=1, start_year=2025, d=30, device_durability_years=30
        )
        self.assertEqual(float(model._get_replacement_cost()), 0)

    def test_invalid_decline_is_rejected(self):
        with self.assertRaises(ValueError):
            storage_service.amazon(drop_rate=1)


if __name__ == "__main__":
    unittest.main()
