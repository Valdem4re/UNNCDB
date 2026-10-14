-- ============================================================
-- Схема БД для датасетов по симметрии кристаллов
-- ============================================================

DROP TABLE IF EXISTS atom_adp CASCADE;
DROP TABLE IF EXISTS atoms CASCADE;
DROP TABLE IF EXISTS structures CASCADE;

-- ------------------------------------------------------------
-- Структуры
-- ------------------------------------------------------------
CREATE TABLE structures (
    id                      BIGSERIAL PRIMARY KEY,
    refcode                 TEXT,
    formula                 TEXT,

    -- Симметрия
    space_group             TEXT,           -- H-M символ как в CIF
    space_group_num         INTEGER,        -- IT номер (1..230), может быть NULL
    point_group             TEXT,           -- 32 точечные группы
    crystal_system          TEXT,           -- 7 сингоний
    is_centrosymmetric      BOOLEAN,        -- TRUE, если есть -1

    -- Ячейка
    a                       DOUBLE PRECISION,
    b                       DOUBLE PRECISION,
    c                       DOUBLE PRECISION,
    alpha                   DOUBLE PRECISION,
    beta                    DOUBLE PRECISION,
    gamma                   DOUBLE PRECISION,
    volume                  DOUBLE PRECISION,
    reduced_volume          DOUBLE PRECISION,  -- volume / Z

    -- Прочее
    z                       INTEGER,
    measurement_temperature DOUBLE PRECISION,
    r_factor_gt             DOUBLE PRECISION,
    wr_factor_gt            DOUBLE PRECISION,
    source_file             TEXT,
    created_at              TIMESTAMPTZ DEFAULT NOW()
);

-- ------------------------------------------------------------
-- Атомы
-- ------------------------------------------------------------
CREATE TABLE atoms (
    structure_id  BIGINT NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    label         TEXT   NOT NULL,
    element       TEXT   NOT NULL,
    x             DOUBLE PRECISION NOT NULL,
    y             DOUBLE PRECISION NOT NULL,
    z             DOUBLE PRECISION NOT NULL,
    u_iso         DOUBLE PRECISION,
    occupancy     DOUBLE PRECISION DEFAULT 1.0,
    PRIMARY KEY (structure_id, label)
);

-- ------------------------------------------------------------
-- Анизотропные ADP
-- ------------------------------------------------------------
CREATE TABLE atom_adp (
    structure_id  BIGINT NOT NULL,
    label         TEXT   NOT NULL,
    u11 DOUBLE PRECISION,
    u22 DOUBLE PRECISION,
    u33 DOUBLE PRECISION,
    u12 DOUBLE PRECISION,
    u13 DOUBLE PRECISION,
    u23 DOUBLE PRECISION,
    PRIMARY KEY (structure_id, label),
    FOREIGN KEY (structure_id, label)
        REFERENCES atoms(structure_id, label) ON DELETE CASCADE
);

-- ------------------------------------------------------------
-- Индексы
-- ------------------------------------------------------------
CREATE INDEX idx_structures_formula          ON structures(formula);
CREATE INDEX idx_structures_space_group      ON structures(space_group);
CREATE INDEX idx_structures_space_group_num  ON structures(space_group_num);
CREATE INDEX idx_structures_centro           ON structures(is_centrosymmetric);
CREATE INDEX idx_structures_volume           ON structures(volume);
CREATE INDEX idx_structures_refcode          ON structures(refcode);
CREATE INDEX idx_atoms_element               ON atoms(element);
CREATE INDEX idx_atoms_structure             ON atoms(structure_id);

-- ------------------------------------------------------------
-- Представление
-- ------------------------------------------------------------
CREATE OR REPLACE VIEW structures_summary AS
SELECT
    s.id,
    s.refcode,
    s.formula,
    s.space_group,
    s.space_group_num,
    s.point_group,
    s.crystal_system,
    s.is_centrosymmetric,
    s.a, s.b, s.c,
    s.alpha, s.beta, s.gamma,
    s.volume,
    s.reduced_volume,
    s.z,
    COUNT(a.label) AS n_atoms
FROM structures s
LEFT JOIN atoms a ON a.structure_id = s.id
GROUP BY s.id;