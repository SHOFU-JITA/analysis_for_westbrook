import requests
import pandas as pd
import time
import os

# 1. 画像から抽出した4地点の緯度・経度リスト
locations = {
    "OL": {"lat": 42.43792, "lon": -72.67160},
    "OS": {"lat": 42.43562, "lon": -72.66720},
    "IL": {"lat": 42.43250, "lon": -72.67651},
    "WB": {"lat": 42.43703, "lon": -72.67615}
}

# 2. 取得したい期間を設定
start_date = "1997-01-01"
end_date = "2025-12-31"

# Open-MeteoのAPIエンドポイント
url = "https://archive-api.open-meteo.com/v1/archive"

# 保存先フォルダ（指定されたパスを使用）
save_dir = r"C:\Users\shoah\codes\Westbrook_P"
# フォルダが存在しない場合は作成する
os.makedirs(save_dir, exist_ok=True)

print("データ取得を開始します...")

# 3. 各地点について順番にAPIリクエストを送信
for name, coords in locations.items():
    print(f"\n[{name}] のデータを取得中...")
    
    # 時間単位（hourly）で気温と降水量を指定
    params = {
        "latitude": coords["lat"],
        "longitude": coords["lon"],
        "start_date": start_date,
        "end_date": end_date,
        "hourly": ["temperature_2m", "precipitation"], # 1時間の気温と降水量
        "timezone": "America/New_York"
    }
    
    response = requests.get(url, params=params)
    
    # データが正常に取得できたか確認
    if response.status_code == 200:
        data = response.json()
        
        # 4. データをPandasのデータフレームに変換
        df = pd.DataFrame({
            "Time": data["hourly"]["time"],
            "Temperature_C": data["hourly"]["temperature_2m"],
            "Precipitation_mm": data["hourly"]["precipitation"]
        })
        
        # 5. CSVファイルとして保存（地点名をファイル名に含める）
        save_path = os.path.join(save_dir, f"{name}_hourly_weather.csv")
        df.to_csv(save_path, index=False, encoding="utf-8-sig")
        
        print(f"成功: CSVファイルの保存が完了しました！\n保存先: {save_path}")
        print("【取得データの一部】")
        print(df.head())
        
    else:
        print(f"エラーが発生しました。地点: {name}, ステータスコード: {response.status_code}")
        
    # 連続でリクエストを送るとサーバーに負荷がかかるため、少し待機時間を設ける
    time.sleep(2)

print("\nすべての地点のデータ取得が完了しました！")