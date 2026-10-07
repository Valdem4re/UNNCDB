-- ============================================================
-- Схема базы данных для хранения CIF-файлов и информации о кристаллических структурах
-- ============================================================

DROP TABLE IF EXISTS atoms CASCADE;
DROP TABLE IF EXISTS structures CASCADE;

-- ------------------------------------------------------------
-- Таблица метаданных структур
-- ------------------------------------------------------------
CREATE TABLE structures (
    id            BIGSERIAL PRIMARY KEY,
    refcode       TEXT UNIQUE,                -- CSD refcode
    formula       TEXT,                       -- Химическая формула
    space_group   TEXT,                       -- Пространственная группа (H-M)
    a             DOUBLE PRECISION,
    b             DOUBLE PRECISION,
    c             DOUBLE PRECISION,
    alpha         DOUBLE PRECISION,
    beta          DOUBLE PRECISION,
    gamma         DOUBLE PRECISION,
    volume        DOUBLE PRECISION,           -- Объём ячейки
    z             INTEGER,                    -- Число формульных единиц
    source_file   TEXT,                       -- Имя исходного CIF-файла
    created_at    TIMESTAMPTZ DEFAULT NOW()
);

-- ------------------------------------------------------------
-- Таблица атомов
-- ------------------------------------------------------------
CREATE TABLE atoms (
    structure_id  BIGINT NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    label         TEXT   NOT NULL,            -- Метка атома (C1, O2, ...)
    element       TEXT   NOT NULL,            -- Химический элемент
    x             DOUBLE PRECISION NOT NULL,
    y             DOUBLE PRECISION NOT NULL,
    z             DOUBLE PRECISION NOT NULL,
    occupancy     DOUBLE PRECISION DEFAULT 1.0,
    u_iso         DOUBLE PRECISION,
    PRIMARY KEY (structure_id, label)
);

-- ------------------------------------------------------------
-- Индексы для быстрого поиска
-- ------------------------------------------------------------
CREATE INDEX idx_structures_formula      ON structures(formula);
CREATE INDEX idx_structures_space_group  ON structures(space_group);
CREATE INDEX idx_structures_volume       ON structures(volume);
CREATE INDEX idx_structures_refcode      ON structures(refcode);
CREATE INDEX idx_atoms_element           ON atoms(element);
CREATE INDEX idx_atoms_structure         ON atoms(structure_id);

-- ------------------------------------------------------------
-- Полезное представление: структуры с числом атомов
-- ------------------------------------------------------------
CREATE OR REPLACE VIEW structures_summary AS
SELECT
    s.id,
    s.refcode,
    s.formula,
    s.space_group,
    s.a, s.b, s.c,
    s.alpha, s.beta, s.gamma,
    s.volume,
    s.z,
    COUNT(a.label) AS n_atoms
FROM structures s
LEFT JOIN atoms a ON a.structure_id = s.id
GROUP BY s.id;

-- ============================================================
-- Готово.
-- ============================================================