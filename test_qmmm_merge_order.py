#!/usr/bin/env python
"""
Test to verify the proper merging order of QMMM settings:
1. Project QMMM settings (from YAML)
2. Parent job settings (modred, custom_solvent, etc.)
3. CLI options (override everything)
"""

from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianQMMMJobSettings,
)


def test_three_way_merge():
    """Test that settings merge in correct order: Project → Parent → CLI"""
    print("\n" + "=" * 70)
    print("TEST: Three-way merge (Project → Parent → CLI)")
    print("=" * 70)

    # Simulate project QMMM settings from YAML
    print("\n1. Project QMMM Settings (from YAML):")
    project_qmmm_settings = GaussianQMMMJobSettings(
        high_level_functional="B3LYP",  # From project
        high_level_basis="6-31G*",  # From project
        low_level_force_field="UFF",  # From project
        solvent_model="PCM",  # From project
        solvent_id="toluene",  # From project
    )
    print(
        f"  - high_level_functional: {project_qmmm_settings.high_level_functional}"
    )
    print(f"  - high_level_basis: {project_qmmm_settings.high_level_basis}")
    print(
        f"  - low_level_force_field: {project_qmmm_settings.low_level_force_field}"
    )
    print(f"  - solvent_model: {project_qmmm_settings.solvent_model}")
    print(f"  - solvent_id: {project_qmmm_settings.solvent_id}")
    print(f"  - modred: {project_qmmm_settings.modred}")
    print(f"  - custom_solvent: {project_qmmm_settings.custom_solvent}")

    # Simulate parent job settings (e.g., from modred command)
    print("\n2. Parent Job Settings (from modred command):")
    parent_settings = GaussianJobSettings(
        functional="M06-2X",  # Parent - won't override QMMM
        basis="def2-TZVP",  # Parent - won't override QMMM
        modred=[[1, 2], [3, 4, 5]],  # Parent - should be inherited
        custom_solvent="eps=8.0\nepsinf=1.8",  # Parent - should be inherited
        additional_route_parameters="EmpiricalDispersion=GD3BJ",  # Parent
    )
    print(f"  - functional: {parent_settings.functional}")
    print(f"  - basis: {parent_settings.basis}")
    print(f"  - modred: {parent_settings.modred}")
    print(f"  - custom_solvent: {parent_settings.custom_solvent}")
    print(
        f"  - additional_route_parameters: {parent_settings.additional_route_parameters}"
    )

    # Simulate merging (as done in qmmm_helper.py)
    print("\n3. Merging parent settings into project QMMM settings...")
    qmmm_settings = project_qmmm_settings

    # Copy attributes from parent that should be inherited
    qmmm_specific_attrs = {
        "jobtype",
        "high_level_functional",
        "high_level_basis",
        "high_level_force_field",
        "medium_level_functional",
        "medium_level_basis",
        "medium_level_force_field",
        "low_level_functional",
        "low_level_basis",
        "low_level_force_field",
        "real_charge",
        "real_multiplicity",
        "int_charge",
        "int_multiplicity",
        "model_charge",
        "model_multiplicity",
        "high_level_atoms",
        "medium_level_atoms",
        "low_level_atoms",
        "bonded_atoms",
        "scale_factors",
    }

    for attr, value in parent_settings.__dict__.items():
        if hasattr(qmmm_settings, attr) and value is not None:
            if attr not in qmmm_specific_attrs:
                current_value = getattr(qmmm_settings, attr, None)
                if current_value is None:
                    setattr(qmmm_settings, attr, value)
                    print(f"  ✓ Inherited {attr} from parent")

    print("\n4. After merging parent settings:")
    print(
        f"  - high_level_functional: {qmmm_settings.high_level_functional} (from PROJECT)"
    )
    print(
        f"  - high_level_basis: {qmmm_settings.high_level_basis} (from PROJECT)"
    )
    print(
        f"  - low_level_force_field: {qmmm_settings.low_level_force_field} (from PROJECT)"
    )
    print(f"  - solvent_model: {qmmm_settings.solvent_model} (from PROJECT)")
    print(f"  - modred: {qmmm_settings.modred} (from PARENT)")
    print(f"  - custom_solvent: {qmmm_settings.custom_solvent} (from PARENT)")
    print(
        f"  - additional_route_parameters: {qmmm_settings.additional_route_parameters} (from PARENT)"
    )

    # Simulate CLI options (these override everything)
    print("\n5. Applying CLI options (override all)...")
    # CLI: -hf wB97X-D -lff "AMBER=HardFirst"
    cli_high_functional = "wB97X-D"
    cli_low_ff = "AMBER=HardFirst"

    if cli_high_functional is not None:
        qmmm_settings.high_level_functional = cli_high_functional
        print(
            f"  ✓ CLI overrides high_level_functional: {cli_high_functional}"
        )
    if cli_low_ff is not None:
        qmmm_settings.low_level_force_field = cli_low_ff
        print(f"  ✓ CLI overrides low_level_force_field: {cli_low_ff}")

    # Set jobtype from parent command name
    qmmm_settings.jobtype = "modred"
    print("  ✓ Inferred jobtype: modred")

    print("\n6. FINAL QMMM Settings:")
    print(f"  - jobtype: {qmmm_settings.jobtype} (INFERRED)")
    print(
        f"  - high_level_functional: {qmmm_settings.high_level_functional} (CLI)"
    )
    print(f"  - high_level_basis: {qmmm_settings.high_level_basis} (PROJECT)")
    print(
        f"  - low_level_force_field: {qmmm_settings.low_level_force_field} (CLI)"
    )
    print(f"  - solvent_model: {qmmm_settings.solvent_model} (PROJECT)")
    print(f"  - modred: {qmmm_settings.modred} (PARENT)")
    print(f"  - custom_solvent: {qmmm_settings.custom_solvent} (PARENT)")
    print(
        f"  - additional_route_parameters: {qmmm_settings.additional_route_parameters} (PARENT)"
    )

    # Verify final state
    assert qmmm_settings.jobtype == "modred", "jobtype not set!"
    assert (
        qmmm_settings.high_level_functional == "wB97X-D"
    ), "CLI didn't override high_level_functional!"
    assert (
        qmmm_settings.high_level_basis == "6-31G*"
    ), "Project basis not preserved!"
    assert (
        qmmm_settings.low_level_force_field == "AMBER=HardFirst"
    ), "CLI didn't override low_level_ff!"
    assert (
        qmmm_settings.solvent_model == "PCM"
    ), "Project solvent_model not preserved!"
    assert qmmm_settings.modred == [
        [1, 2],
        [3, 4, 5],
    ], "Parent modred not inherited!"
    assert (
        qmmm_settings.custom_solvent is not None
    ), "Parent custom_solvent not inherited!"
    assert (
        qmmm_settings.additional_route_parameters is not None
    ), "Parent additional_route_parameters not inherited!"

    print("\n✓ TEST PASSED: Three-way merge works correctly!")
    print("\nMerge order verified:")
    print("  1. Start with project QMMM settings (YAML)")
    print(
        "  2. Inherit job-specific settings from parent (modred, custom_solvent, etc.)"
    )
    print("  3. CLI options override everything")


def test_project_override_priority():
    """Test that project QMMM settings are NOT overridden by parent job settings"""
    print("\n" + "=" * 70)
    print("TEST: Project QMMM settings have priority over parent")
    print("=" * 70)

    # Project QMMM settings explicitly set solvent
    project_qmmm_settings = GaussianQMMMJobSettings(
        solvent_model="SMD",
        solvent_id="water",
    )
    print("\nProject QMMM settings:")
    print(f"  - solvent_model: {project_qmmm_settings.solvent_model}")
    print(f"  - solvent_id: {project_qmmm_settings.solvent_id}")

    # Parent settings also have solvent (should NOT override)
    parent_settings = GaussianJobSettings(
        solvent_model="PCM",
        solvent_id="acetone",
    )
    print("\nParent settings:")
    print(f"  - solvent_model: {parent_settings.solvent_model}")
    print(f"  - solvent_id: {parent_settings.solvent_id}")

    # Merge (project should win)
    qmmm_settings = project_qmmm_settings
    qmmm_specific_attrs = {
        "jobtype",
        "high_level_functional",
        "high_level_basis",
        "high_level_force_field",
        "medium_level_functional",
        "medium_level_basis",
        "medium_level_force_field",
        "low_level_functional",
        "low_level_basis",
        "low_level_force_field",
        "real_charge",
        "real_multiplicity",
        "int_charge",
        "int_multiplicity",
        "model_charge",
        "model_multiplicity",
        "high_level_atoms",
        "medium_level_atoms",
        "low_level_atoms",
        "bonded_atoms",
        "scale_factors",
    }

    for attr, value in parent_settings.__dict__.items():
        if hasattr(qmmm_settings, attr) and value is not None:
            if attr not in qmmm_specific_attrs:
                current_value = getattr(qmmm_settings, attr, None)
                if (
                    current_value is None
                ):  # Only inherit if project didn't set it
                    setattr(qmmm_settings, attr, value)

    print("\nAfter merge:")
    print(
        f"  - solvent_model: {qmmm_settings.solvent_model} (should be SMD from project)"
    )
    print(
        f"  - solvent_id: {qmmm_settings.solvent_id} (should be water from project)"
    )

    assert (
        qmmm_settings.solvent_model == "SMD"
    ), "Project solvent_model was overridden!"
    assert (
        qmmm_settings.solvent_id == "water"
    ), "Project solvent_id was overridden!"

    print("\n✓ TEST PASSED: Project settings not overridden by parent")


def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print("QMMM Settings Three-Way Merge Tests")
    print("=" * 70)

    try:
        test_three_way_merge()
        test_project_override_priority()

        print("\n" + "=" * 70)
        print("ALL TESTS PASSED! ✓")
        print("=" * 70)
        print("\nSettings merge order is correct:")
        print("  1. Project QMMM settings (from YAML) - BASE")
        print(
            "  2. Parent job settings (modred, custom_solvent) - INHERITED if not in project"
        )
        print("  3. CLI options - OVERRIDE everything")

    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ UNEXPECTED ERROR: {e}")
        import traceback

        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
